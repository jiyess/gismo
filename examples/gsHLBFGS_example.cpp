/* @file -

@brief - A reference implementation of the following paper.
	If you make use of the code or the idea/algorithm in your work, please cite our paper
	Ji, Y., Yu, Y. Y., Wang, M. Y., & Zhu, C. G. (2021).
	Constructing high-quality planar NURBS parameterization for
	isogeometric analysis by adjustment control points and weights.
	Journal of Computational and Applied Mathematics, 396, 113615.
	(https://www.sciencedirect.com/science/article/pii/S0377042721002375)

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): Ye Ji (jiye@mail.dlut.edu.cn)
*/

#include <gismo.h>
#include <gsHLBFGS/gsHLBFGS.h>
//#include <gsModeling/gsPatchGenerator.h>

using namespace gismo;

template<typename T>
gsVector<T> convert_mp_to_gsFreeVec(const gsMultiPatch<T> &mp,
                                    const gsDofMapper &mapper) {

    gsVector<T> freeVec(mapper.freeSize());

    for (index_t iptch = 0; iptch != mp.nPatches(); iptch++) {
        for (index_t idim = 0; idim != mp.targetDim(); idim++) {
            for (index_t idof = 0; idof != mapper.patchSize(iptch, idim); idof++)
                // if it is possible to just loop over the free index
                // since this function is called very often during optimization
                if (mapper.is_free(idof, iptch, idim)) {
                    index_t idx = mapper.index(idof, iptch, idim);
                    freeVec(idx) = mp.patch(iptch).coefs()(idof, idim);
                }
        }
    }
    return freeVec;
}

template<typename T>
void convert_gsFreeVec_to_mp(const gsVector<T> &gsFreeVec,
                             const gsDofMapper &mapper,
                             gsMultiPatch<T> &mp) {
    for (index_t iptch = 0; iptch != mp.nPatches(); iptch++) {
        for (index_t idim = 0; idim != mp.targetDim(); idim++) {
            for (index_t idof = 0; idof != mapper.patchSize(iptch, idim); idof++)
                // if it is possible to just loop over the free index
                // since this function is called very often during optimization
                if (mapper.is_free(idof, iptch, idim)) {
                    index_t idx = mapper.index(idof, iptch, idim);
                    mp.patch(iptch).coefs()(idof, idim) = gsFreeVec(idx);
                }
        }
    }
}

template<short_t d, typename T>
class gsObjInterFreeFunc : public gsOptProblem<T> {
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space space;
    typedef gsExprAssembler<>::solution solution;

    using Base = gsOptProblem<T>;
public:
    gsObjInterFreeFunc(const gsMultiPatch<T> &patches,
                       gsDofMapper mapper,
                       const T eps
    )
            :
            m_mp(patches),
            m_mapper(std::move(mapper)),
            m_mb(m_mp),
            m_eps(eps) {
        m_assembler.setIntegrationElements(m_mb);
        space u = m_assembler.getSpace(m_mb, d); // 1D space!!
        m_evaluator = gsExprEvaluator<T>(m_assembler);
    }

    T evalObj(const gsAsConstVector<T> &u) const final {
        gsMultiPatch<T> mp = m_mp;
        convert_gsFreeVec_to_mp<T>(u, m_mapper, mp);
        geometryMap G = m_evaluator.getMap(mp);

        auto EfoldElimination = (m_eps - jac(G).det()).ppartval();

        T F = m_evaluator.integral(EfoldElimination);
        return F;
    }

    void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const final {
        gsMultiPatch<T> mp = m_mp;
        convert_gsFreeVec_to_mp<T>(u, m_mapper, mp);
        geometryMap G = m_assembler.getMap(mp);

        space space = m_assembler.getSpace(m_mb, d); // 1D space!!
        space.setupMapper(m_mapper);

//      |J|' w.r.t. physical coordinates x and y
        auto derJacDet = jac(space) % jac(G).tr().adj();

        // TODO:下面这行先这么写，回头修改掉(出于效率考虑)，引入新的expression？
        auto indicator = -heaviside(m_eps - jac(G).det());
        auto E_der = indicator * derJacDet;

        m_assembler.initSystem();
        m_assembler.assemble(E_der);
        result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()), m_assembler.rhs().rows());
    }

protected:
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsDofMapper m_mapper;

    mutable gsExprEvaluator<T> m_evaluator;
    mutable gsExprAssembler<T> m_assembler;

    T m_eps; // need to handle later, set m_eps = 0.05*S
};

template<short_t d, typename T>
class gsObjQualityImprovePt : public gsOptProblem<T> {
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space space;
    typedef gsExprAssembler<>::solution solution;

    using Base = gsOptProblem<T>;
public:
    gsObjQualityImprovePt(const gsMultiPatch<T> &patches,
                          gsDofMapper mapper,
                          const T lambda1,
                          const T lambda2)
            :
            m_mp(patches),
            m_mapper(std::move(mapper)),
            m_mb(m_mp),
            m_lambda1(lambda1),
            m_lambda2(lambda2) {
        m_assembler.setIntegrationElements(m_mb);
        space u = m_assembler.getSpace(m_mb, d); // 1D space!!
        m_evaluator = gsExprEvaluator<T>(m_assembler);
    }

    T evalObj(const gsAsConstVector<T> &u) const final {
        gsMultiPatch<T> mp = m_mp;
        convert_gsFreeVec_to_mp<T>(u, m_mapper, mp);

        geometryMap G = m_evaluator.getMap(mp);

        gsConstantFunction<T> lambda1(m_lambda1, d);
        gsConstantFunction<T> lambda2(m_lambda2, d);
        auto lam1 = m_evaluator.getVariable(lambda1);
        auto lam2 = m_evaluator.getVariable(lambda2);

        gsConstantFunction<T> penaltyFactor(1e20, d);
        auto penaF = m_evaluator.getVariable(penaltyFactor);

        auto Ewinslow = jac(G).sqNorm() / jac(G).det();
        auto Euniform = pow(jac(G).det(), 2);
        auto EimprovePtReal = lam1 * Ewinslow + lam2 * Euniform;

        // TODO: need to compute two times now? how to improve it?
        auto indicator = heaviside(jac(G).det());
        auto indicatorNot = heaviside(-jac(G).det());
        auto EimprovePt = indicator * EimprovePtReal + indicatorNot * penaF;

        T F = m_evaluator.integral(EimprovePt);
        return F;
    }

    void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const final {
        gsMultiPatch<T> mp = m_mp;
        convert_gsFreeVec_to_mp<T>(u, m_mapper, mp);
        geometryMap G = m_assembler.getMap(mp);

        space space = m_assembler.getSpace(m_mb, d); // 1D space!!
        space.setupMapper(m_mapper);

//      |J|' w.r.t. physical coordinates x and y
        auto derJacDet = jac(space) % jac(G).tr().adj();

        // TODO: are there a mass of redundant computation. If yes, how to avoid it?
        auto Ewinslow = jac(G).sqNorm() / jac(G).det();
        auto derEwinslow = 2.0 / jac(G).det() * (jac(space) % jac(G)) - Ewinslow.val() / jac(G).det() * derJacDet;
        auto derEuniform = 2.0 * jac(G).det() * derJacDet;

        auto derEimprovePt = m_lambda1 * derEwinslow + m_lambda2 * derEuniform;

        m_assembler.initSystem();
        m_assembler.assemble(derEimprovePt);
        result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()), m_assembler.rhs().rows());
    }

protected:
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsDofMapper m_mapper;

    mutable gsExprEvaluator<T> m_evaluator;
    mutable gsExprAssembler<T> m_assembler;

    T m_lambda1, m_lambda2;
};

// It is a barrier function-based method. The main step ff
template<short d, typename T>
class gsBarrierMethod : public gsPatchGenerator<T> {
//    typedef typename gsExprEvaluator<T>::geometryMap geometryMap;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space space;
    typedef gsExprAssembler<>::solution solution;
public:

    gsBarrierMethod(const gsMultiPatch<T> &bRep,
                    const index_t &method, const bool &plot_init)
            :
            m_bRep(bRep),
            m_method(method),
            m_plot_init(plot_init) {
        GISMO_ASSERT(d == 2, "This method is only available for d==2");
        // TODO: assign m_mp by using initialization

        //sanityCheckInput

        initialization();
        m_mb = gsMultiBasis<T>(m_mp);
        makeMapper();

        m_area = computeArea();
        m_eps = 0.05 * m_area; // ATTENTION!!
    }

public:

//    const gsGeometry<T> &compute() {GISMO_NO_IMPLEMENTATION; }
//    void enablePlot() { m_plot = true; }
//    void disablePlot() { m_plot = false; }

    gsDofMapper mapper() { return m_mapper; }

    void initialization() {
        // TO DO: different initialization methods:
        // 1. discrete Coons 2. Smoothness energy 3. Spring model etc.
        switch (m_method) {
            case 1: {
                // Spring model method
                gsInfo << "Using spring patch construction.\n";
                gsSpringPatch<T> spring(m_bRep);
                gsInfo << "Created a " << spring.compute() << "\n";
                //if (save) gsWrite(spring.result(), "result_patch");
                m_mp.addPatch(spring.result());
                if (m_plot_init)
                    gsWriteParaview(m_mp, "mp_init_spring", 1000, true, true);
                break;
            }
            case 2: {
                //// Cross Approximation patch method
                //// Question: only works for 3D surfaces? i.e., for mapping: x: \mathbb{R}^2 --> \mathbb{R}^3
                //// I am not sure, ask Hugo or Matthias later.

                //gsInfo << "Using cross approximation construction.\n";
                //gsCrossApPatch<T> cross(m_bRep);
                //gsDebug << "xxxxxxxxxxxxxxxxxxx" << "\n";
                //gsInfo << "Created a " << cross.compute() << "\n";
                ////if (save) gsWrite(spring.result(), "result_patch");
                //m_mp.addPatch(cross.result());
                //if (m_plot_init)
                //	gsWriteParaview(m_mp, "mp_init_cross", 1000, true, true);
                //break;
            }
            case 3: {
                // consturt a parameterization with the inner control points all equal to (0, 0)

                gsInfo << "Set all the inner control points to origin of coordinates.\n";
                gsCoonsPatch<T> coons(m_bRep);
                coons.compute();
                m_mp.addPatch(coons.result());

                gsVector<T> initialZeros(m_mapper.freeSize());
                convert_gsFreeVec_to_mp(initialZeros, m_mapper, m_mp);
                gsInfo << "Created a same point Patch." << "\n";

            }
            case 4: {
                // Smoothness energy method
                // TODO: need to implement later...
                // However, the results seems the same as Spring model method?

                break;
            }
            case 0:
            default:
                gsInfo << "Using Coons' patch construction.\n";
                gsCoonsPatch<T> coons(m_bRep);
                gsInfo << "Created a " << coons.compute() << "\n";
                //if (save) gsWrite(coons.result(), "result_patch");
                m_mp.addPatch(coons.result());
                if (m_plot_init)
                    gsWriteParaview(m_mp, "mp_init_coons", 1000, true, true);
                break;
        }
    }

    void makeMapper() {
        // Now, we set all the inner control points as optimization variables
        // It is possible to set only a part of them as optimization variables later

        ////! [Make mapper for the design DoFs]
        m_mapper.init(gsMultiBasis<>(m_mp), m_mp.targetDim());
        //gsDofMapper m_mapper(gsMultiBasis<>(m_mp), m_mp.targetDim());
        // 1. Mark the vertical displacements of the control points (except the corners)
        //    as design variables
        //      1.1. Eliminate boundary control points

        gsMatrix<index_t> idx;
        for (index_t iptch = 0; iptch != m_mp.nPatches(); iptch++) {
            idx = m_mp.basis(iptch).allBoundary(); // if it needs to compute all basis or not?
            // if YES, is there any more efficient way?
            for (index_t idim = 0; idim != m_mp.targetDim(); idim++) {
                for (index_t idof = 0; idof != idx.size(); idof++) {
                    m_mapper.eliminateDof(idx(idof), iptch, idim);
                }
            }
        }
        m_mapper.finalize();

        gsInfo << "#Numb of free  variables is " << m_mapper.freeSize() << "\n";
        gsInfo << "#Numb of fixed variables is " << m_mapper.boundarySize() << "\n";
        gsInfo << "#Numb of total variables is " << m_mapper.size() << "\n";
    }

    T computeArea() {
        // compute the area of the computational domain by B-Rep
        // using the following Green's formulation:
        // S = \int_{\Omega} 1 d \Omega
        //   = \oint_{\partial \Omega} x(t) y'(t) dt
        //   = \sum_{1}^4 \int_0^1 x_i(t) y'_i(t) dt

        // Here, one must take care of the orientation of the boundary curves
        // I found there exist some files get negative values of area
        // We should make this part more robust!! make it indenpent of boundary orientation
        // make some pre-check? or just get the absolute value of the result?

        // Or if the opposite (like NS and WE) boundary curves always with the same direction?

        gsMultiPatch<T> mp = m_mp;
        m_evaluator.setIntegrationElements(m_mb);
        geometryMap G = m_evaluator.getMap(mp);

        // TODO: 由于初始化通常折叠的问题，这种方法的计算结果非常不稳定，先凑合用; 利用文章中沿着边界的积分来计算
        T area = m_evaluator.integral(jac(G).det());

        return area;
    }

    const gsGeometry<T> &compute() final {
        /////////////////////////// STEP 2: Foldovers Elimination Step ///////////////////////////////
        gsVector<T> initU = convert_mp_to_gsFreeVec<T>(m_mp, m_mapper);

        gsObjInterFreeFunc<2, real_t> obj(m_mp, m_mapper, m_eps);

        gsHLBFGS<real_t> optimizer(&obj);
        // check inital guess
        gsStopwatch stopwatch;
        optimizer.solve(initU);

        ///////////////////// STEP 3: improving the quality of parameterization ///////////////////////////
        //// Another optimization problem
        gsObjQualityImprovePt<2, real_t> objImprove(m_mp, m_mapper, 1.0, 1.0 / std::pow(m_area, 2));

        gsHLBFGS<real_t> optimizer2(&objImprove);
        optimizer2.solve(optimizer.currentDesign());
        gsInfo << "running time : " << stopwatch.stop() << "\n";

        convert_gsFreeVec_to_mp<T>(optimizer2.currentDesign(), m_mapper, m_mp);
        //--------------------------------------------------------------------------------------------------

        outputResult();
        return m_mp.patch(0);
    }

    void outputResult() const {

        m_evaluator.setIntegrationElements(m_mb);
        geometryMap G = m_evaluator.getMap(m_mp);
        GISMO_ENSURE(m_mp.nPatches() == 1, "Does not yet work for multipatch, "
                                           "but multipatch has " << m_mp.nPatches() << " patches");

        // resulting mesh
        index_t numSamples(1000);

        bool plot_mesh(false);
        bool plot_net(false);
        gsWriteParaview(m_mp, "mp", numSamples, plot_mesh, plot_net);

        // Scaled Jacobian
        m_evaluator.options().setInt("plot.npts", numSamples);
        auto metric_ScaledJacobian = jac(G).det() / (jac(G)[0].norm() * jac(G)[1].norm());
        m_evaluator.writeParaview(metric_ScaledJacobian, G, "metric_ScaledJacobian");

        // Frobenius condition number
        auto metric_FrobCondNum = jac(G).sqNorm() / meas(G);
        m_evaluator.writeParaview(metric_FrobCondNum, G, "metric_FrobCondNum");
    }

private:

    gsDofMapper m_mapper;
    mutable gsExprEvaluator<T> m_evaluator;
    mutable gsExprAssembler<T> m_assembler;

    mutable gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    const gsMultiPatch<T> m_bRep;

    T m_area; // area of computational domain
    T m_eps; // parameter need for foldover elimination step

    index_t m_method;
    bool m_plot_init;
};

int main(int argc, char *argv[]) {

    //////////////////// STEP 1: read a boundary represention (B-Rep) file /////////////////
    bool save = false;
    index_t method = 0;
    real_t tol = 1e-10;
    bool plot_init = false;

    // Load XML file containing the boundary represention
    // TODO: give some different cases as defalut benchmarks,
    //       otherwise deal with the input filename
    std::string filename_input("breps/2D/duck_boundary.xml");
//    std::string filename_input("breps/2D/puzzle4_bdry.xml");
    std::string filename_output("results");

    // Read input from command line arguments
    //! [Parse command line]
    gsCmdLine cmd("Hi, give me a file (eg: .xml) containing boundary representation (B-Rep)"
                  "and I will try to parameterize it!");

    cmd.addPlainString("input", "Name of the input file containing boundary data", filename_input);
    cmd.addString("o", "output", "Name of the output file", filename_output);
    cmd.addInt("m", "method", "Method: 0 Coons' patch (default), 1 Spring patch, 2: Cross-Ap. patch", method);
    cmd.addReal("t", "tolerance", "Tolerance for identifing patch interfaces", tol);
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("plotInit", "Plot resulting initaialization for Paraview", plot_init);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }
    //! [Parse command line]

    // Load XML file
    gsMultiPatch<real_t> bRep;
    gsReadFile<>(filename_input, bRep);
    GISMO_ENSURE(!bRep.empty(), "The gsMultiPatch is empty - maybe file is missing or corrupt.");

    // TODO: At present， we need the input B-Rep is ordered by (W->E->N->S),
    //  make it more robust and available to whatever the order of B-Rep curves

    /////////////////////////////////////// STEP 2: eliminate foldovers ////////////////////////////////
    gsBarrierMethod<2, real_t> opt(bRep, method, plot_init);
    opt.compute();

    // TODO: scale the computational domain for better numerical stability, using the area?

    ///////////////////////////// STEP 4: output the resulting parameterization ////////////////////////
    // writing the resulting parameterization to a G+Smo .xml file
    // filename_output is a string. The extention .xml is added automatically
    // TODO: other formats? To make it easy to visulization

    ////////////////////////////////////// STEP 5: VISUALIZATION ///////////////////////////////////////
    // TODO: visualization, in MATLAB first? then ?
    // GNUPLOT, looks good; a MATLAB plot style

    return EXIT_SUCCESS;
}

