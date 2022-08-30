/** @file gsBarrierPatch.h

    @brief XXXXXXXX

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): XXXXXXXX
*/

#pragma once

#include<gsModeling/gsPatchGenerator.h>
#include<gsAssembler/gsExprAssembler.h>
#include<gsAssembler/gsExprEvaluator.h>
#include<gsOptimizer/gsOptProblem.h>
#include<gsOptimizer/gsOptimizer.h>

namespace gismo {


/**
\brief XXXXXXXX

\tparam XXXXXXXX

\ingroup Modeling
*/
template<short_t d, typename T>
class gsBarrierPatch : public gsPatchGenerator<T> {
public:
    typedef gsPatchGenerator<T> Base;

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space space;
    typedef gsExprAssembler<>::solution solution;

public:

    /** 
       \brief XXXXXXXX

       XXXXXXXX

       \param XXXXXXXX
    */
    gsBarrierPatch(const gsMultiPatch<T> &boundary)
    :
    Base(boundary)
    {

    }

    gsBarrierPatch(const gsMultiPatch<T> &bRep,
                    const index_t &method, const bool &plot_init);


public:

    // Look at gsPatchGenerator
    const gsGeometry<T> &compute();

protected:

    using Base::m_boundary;

    using Base::m_result;

}; // gsBarrierPatch

template<short_t d, typename T>
class gsObjQualityImprovePt : public gsOptProblem<T>
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space space;
    typedef gsExprAssembler<>::solution solution;

    using Base = gsOptProblem<T>;
public:
    gsObjQualityImprovePt(const gsMultiPatch<T> &patches,
                          const gsDofMapper &mapper,
                          const T eps);

    T evalObj(const gsAsConstVector<T> &u) const;

    void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const {
        gsMultiPatch<T> mp = m_mp;
        this->convert_design_to_mp(u, mp);
        geometryMap G = m_assembler.getMap(mp);

        space space = m_assembler.getSpace(m_mb, d); // 1D space!!
        space.setupMapper(m_mapper);

//        gsConstantFunction<T> lambda1(m_lambda1, d);
//        gsConstantFunction<T> lambda2(m_lambda2, d);
//        auto lam1 = m_assembler.getCoeff(lambda1);
//        auto lam2 = m_assembler.getCoeff(lambda2);

//      |J|' w.r.t. physical coordinates x and y
        auto derJacDet = jac(space) % jac(G).tr().adj();

        // TODO:下面这行先这么写，回头修改掉(出于效率考虑)，引入新的expression？
//        auto indicator = (m_eps-jac(G).det()).ppartval()/(1e-12-(m_eps-jac(G).det()).ppartval());
//        auto E_der = indicator * derJacDet;

        auto Ewinslow = (jac(G) % jac(G)) / jac(G).det();
        auto derEwinslow = 2.0 / jac(G).det() * (jac(space) % jac(G)) - Ewinslow.val() / jac(G).det() * derJacDet;
        auto derEuniform = 2.0 * jac(G).det() * derJacDet;

        gsConstantFunction<T> lambda1(m_lambda1, d);
        gsConstantFunction<T> lambda2(m_lambda2, d);
        auto lam1 = m_evaluator.getVariable(lambda1);
        auto lam2 = m_evaluator.getVariable(lambda2);

        gsMatrix<T> PF(1, 1);
        PF << 1e20;
        gsConstantFunction<T> penaltyFactor(PF, d);
        auto penaF = m_evaluator.getVariable(penaltyFactor);

//        gsDebugVar(penaF);

        auto derEimprovePt = m_lambda1 * derEwinslow + m_lambda2 * derEuniform;

        // 这里加标量, 加不上去
//        auto indicator = (-jac(G).det()).ppartval()/(1e-20-(-(jac(G).det()).ppartval()));
//        auto indicator2 = -jac(G).det().ppartval()/(1e-20-((jac(G).det()).ppartval()));
//        auto derEimprovePt = derEimprovePtReal * indicator2  ;

//        auto indicator =  heaviside(jac(G).det());
//        auto indicatorNot = heaviside(-jac(G).det());
//        auto derEimprovePt = derEimprovePtReal*indicator + ;

//        auto derEimprovePt2 = penaF * indicator;
//        auto derEimprovePt = derEimprovePt1+derEimprovePt2;
//
//        gsDebugVar(derEimprovePt.rows());
//        gsDebugVar(derEimprovePt.cols());

//        gsVector<T> pt(2);
//        pt.setConstant(0.5);
//        gsDebugVar(m_evaluator.eval(indicator, pt));
//        gsDebugVar(m_evaluator.eval(indicator2, pt));
//        gsDebugVar(m_evaluator.eval(derEimprovePt, pt));
//        gsDebugVar(m_evaluator.eval(test, pt));
//        gsDebugVar(m_evaluator.eval(E_der, pt));

//     现在看来可能是这里有问题？？跟Hugo讨论一下！
        m_assembler.initSystem();
        m_assembler.assemble(derEimprovePt);
//        gsDebugVar(m_assembler.matrix());
        result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()), m_assembler.rhs().rows());

//        // CHECK INVERSES OF jac(G)
//        auto JTG = (jac(space) % jac(G));//.tr(); // transpose to switch to rowSpace
//        auto JiG = (jac(space) % jac(G).inv());//.tr(); // transpose to switch to rowSpace
//        auto J_prime = jac(G).det() * JiG;
//        auto Ewinslow = (jac(G).tr() * jac(G)).trace() / jac(G).det();
//        auto Ewinslow_der1 = (2.0 / jac(G).det()) * JTG;//
//        auto Ewinslow_der2 = -(Ewinslow.val() / jac(G).det()) * J_prime;
//        auto Ewinslow_der = Ewinslow_der1 + Ewinslow_der2;
//        auto Euniform_der = 2.0 / pow(S, 2) * (jac(G).det().val() - S.val()) * J_prime;
//
//        // To evaluate expressions at the point 'pt', use the following
////        gsVector<T> pt(2);
////        pt.setConstant(0.5);
////        gsDebugVar(m_evaluator.eval(G, pt));
////        gsDebugVar(m_evaluator.eval(jac(space), pt));
////        gsDebugVar(m_evaluator.eval(jac(G), pt));
////        gsDebugVar(m_evaluator.eval(JTG, pt));
////        gsDebugVar(JTG.Space);
////        gsDebugVar(m_evaluator.eval(JiG, pt));
////        gsDebugVar(JiG.Space);
////
////        gsDebugVar(m_evaluator.eval(J_prime, pt));
////
////        gsDebugVar(m_evaluator.eval(Ewinslow, pt));
////        gsDebugVar(m_evaluator.eval(Ewinslow_der1, pt));
////        gsDebugVar(m_evaluator.eval(Ewinslow_der2, pt));
////        gsDebugVar(m_evaluator.eval(Ewinslow_der, pt));
////        gsDebugVar(m_evaluator.eval(Euniform_der, pt));
//
//        m_assembler.initSystem();
//        m_assembler.assemble(m_lambda1 * Ewinslow_der + m_lambda2 * Euniform_der);
//        result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()), m_assembler.rhs().rows());
    }

protected:
    void convert_design_to_mp(const gsAsConstVector<T> &u, gsMultiPatch<T> &mp) const {
        // Make the geometry
        index_t idx;
        for (size_t i = 0; i != mp.nPatches(); i++) {
            for (size_t c = 0; c != mp.targetDim(); c++) {
                for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)

                    // if it is possible to just loop over the free index
                    // since this function is called very often during optimization
                    if (m_mapper.is_free(k, i, c)) {
                        idx = m_mapper.index(k, i, c);
                        mp.patch(i).coefs()(k, c) = u[idx];
                    }
            }
        }
    }

    gsVector<T> convert_mp_to_gsvec(const gsMultiPatch<T> &mp) const {
        // TODO: Now, it just set all free variables into the vector u
        // I will integrate it with the initialization step

        gsVector<T> currentDesign(m_mapper.freeSize());

        index_t idx;
        for (size_t i = 0; i != mp.nPatches(); i++) {
            for (size_t c = 0; c != mp.targetDim(); c++) {
                for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)
                    // if it is possible to just loop over the free index
                    // since this function is called very often during optimization
                    if (m_mapper.is_free(k, i, c)) {
                        idx = m_mapper.index(k, i, c);
                        currentDesign(idx) = mp.patch(i).coefs()(k, c);
                    }
            }
        }

        return currentDesign;
    }

protected:
    gsMultiPatch<T> m_mp;
    gsDofMapper m_mapper;
    index_t m_patchID;

    gsMultiBasis<T> m_mb;

    mutable gsExprEvaluator<T> m_evaluator;
    mutable gsExprAssembler<T> m_assembler;

    gsOptionList m_options;
    T m_lambda1, m_lambda2;


    /// Number of design variables
    using Base::m_numDesignVars;

    /// Number of constraints
    using Base::m_numConstraints;

    /// Number of nonzero entries in the Constraint Jacobian
    using Base::m_numConJacNonZero;

    /// Lower bounds for the design variables
    using Base::m_desLowerBounds;

    /// Upper bounds for the design variables
    using Base::m_desUpperBounds;

    /// Lower bounds for the constraints
    using Base::m_conLowerBounds;

    /// Upper bounds for the constraints
    using Base::m_conUpperBounds;

    /// Constraint Jacobian non-zero entries rows
    using Base::m_conJacRows;

    /// Constraint Jacobian non-zero entries columns
    using Base::m_conJacCols;

    /// Current design variables (and starting point )
    using Base::m_curDesign;

    mutable gsMatrix<T> m_allBasisVal;
    mutable gsMatrix<T> m_allBasis1stDervs;
    mutable gsMatrix<T> m_allBasis2ndDervs;
    mutable gsMatrix<T> m_gaussWts;
    mutable gsMatrix<index_t> m_gaussIdx;

    T m_eps; // need to handle later, set m_eps = 0.05*S
};

template<short_t d, typename T>
class gsObjInterFreeFunc : public gsOptProblem<T> {
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space space;
    typedef gsExprAssembler<>::solution solution;

    using Base = gsOptProblem<T>;
public:
    gsObjInterFreeFunc(const gsMultiPatch<T> &patches,
                       const gsDofMapper &mapper,
                       const T eps
    )
            :
            m_mp(patches),
            m_mapper(mapper),
            m_mb(m_mp),
            m_patchID(0),
            m_eps(eps) {

        // What does the following part mean? seems useless...

        // parameters related with design variables
        // Number of design variables
        m_numDesignVars = m_mapper.freeSize();

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        // parameters related with constraints
        // There is no constraint in our method
        m_numConstraints = 0;

        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);
        //m_conLowerBounds[0] = m_conUpperBounds[0] = 0;

        // parameters related with nonzero entries in the Constraint Jacobian
        // this->currentDesign_into(mp, m_curDesign);
        m_numConJacNonZero = m_numDesignVars;
        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);
        m_curDesign = convert_mp_to_gsvec(m_mp);

        m_assembler.setIntegrationElements(m_mb);
        space u = m_assembler.getSpace(m_mb, d); // 1D space!!
        m_evaluator = gsExprEvaluator<T>(m_assembler);

        m_lambda1 = 1.0;
        m_lambda2 = 0.0;
    }

    T evalObj(const gsAsConstVector<T> &u) const {
        T F = 0;
        gsMultiPatch<T> mp = m_mp;
        this->convert_design_to_mp(u, mp);

        geometryMap G = m_evaluator.getMap(mp);

        //ReLU, Rectified Linear Unit, how to implement?

//        auto jacGdet = jac(G).det();
//        auto EfoldElimination = jacGdet;//.val();
        auto EfoldElimination = (m_eps - jac(G).det()).ppartval();//.val();

//        gsDebugVar(m_eps);

        // 这里需要一个指示函数，目前只对不为0的部分做测试
//        gsConstantFunction<T> lambda1(-0.5, d);
//        auto lam1 = m_evaluator.getVariable(lambda1);
//        auto EfoldElimination = jac(G).det();
//        gsVector<T> pt(2);
//        pt.setConstant(0.3);
//        gsDebugVar(m_evaluator.eval(EfoldElimination, pt));
//        gsDebugVar(m_evaluator.eval(m_eps-jac(G).det(), pt));

//        auto EfoldElimination = jac(G).det() * jac(G).det(); // for test purpose
        F = m_evaluator.integral(EfoldElimination);
//          F=0;
        // Here, what is the meaning?
//        gsConstantFunction<T> lambda1(m_lambda1, d);
//        gsConstantFunction<T> lambda2(m_lambda2, d);
//        auto lam1 = m_evaluator.getVariable(lambda1);
//        auto lam2 = m_evaluator.getVariable(lambda2);

        // Q: here, need some extra calculation?
//        auto Ewinslow = (jac(G).tr() * jac(G)).trace() / jac(G).det();
//        auto Euniform = gismo::expr::pow(jac(G).det(), 2);// or pow(jac(G).det() / S - 1,2) --> what is S?

        // if (d==2)
//        F = m_evaluator.integral(lam1 * Ewinslow.val() + lam2 * Euniform);

        // else if (d==3)

        return F;
    }

    void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const {
        gsMultiPatch<T> mp = m_mp;
        this->convert_design_to_mp(u, mp);
        geometryMap G = m_assembler.getMap(mp);

//        gsConstantFunction<T> lambda1(m_lambda1, d);
//        gsConstantFunction<T> lambda2(m_lambda2, d);
//        auto lam1 = m_assembler.getCoeff(lambda1);
//        auto lam2 = m_assembler.getCoeff(lambda2);

        space space = m_assembler.getSpace(m_mb, d); // 1D space!!
        space.setupMapper(m_mapper);

// 问题应该就在这!!
//        gsMultiBasis<> basis(m_mb);
//        space space = m_assembler.getSpace(basis, 1); // 1D space!!
//        space.setupMapper(m_mapper);

        m_evaluator = gsExprEvaluator<T>(m_assembler);
        T area = m_evaluator.integral(meas(G));
//
//        gsVector<T> pt(2);
//        pt.setConstant(0.5);
//        gsDebugVar(m_evaluator.eval(meas(G), pt));
        gsDebugVar(area);

//        gsConstantFunction<T> area_fun(area, d);
//        auto S = m_assembler.getCoeff(area_fun);

//        gsVector<T> pt(2);
////        pt.setConstant(0.3);
//        pt(0) = 0.177459666924148;
//        pt(1) = 0.141967733539319;
//        gsDebugVar(pt);

        // what is jac(space)?? we need the derivatives of basis functions, not space
//        auto basisDers = jac(space);
//        auto adjJac = jac(G).tr().adj();
//        gsDebugVar(m_evaluator.eval(basisDers, pt));
//        gsDebugVar(m_evaluator.eval(adjJac, pt));
//        auto jacG = jac(G);
//        gsDebugVar(m_evaluator.eval(jacG, pt));
//        auto invJac = jac(G).inv();
//        gsDebugVar(m_evaluator.eval(invJac, pt));

//      |J|' w.r.t. physical coordinates x and y
        auto derJacDet = jac(space) % jac(G).tr().adj();
//        gsDebugVar(m_evaluator.eval(basisDersTest, pt));
//        gsDebugVar(jac(space).rows());
//        gsDebugVar(jac(space).cols());

//      这个有用！！！！
//        auto adjTest = jac(G).adj();
//        gsDebugVar(m_evaluator.eval(adjTest, pt));

//        auto jacTest = jac(G).inv()*jac(G).det();
//        gsDebugVar(m_evaluator.eval(jacTest,pt));

        // TODO:下面这行先这么写，回头修改掉(出于效率考虑)，引入新的expression？
//        auto indicator = (m_eps-jac(G).det()).ppartval()/(1e-12-(m_eps-jac(G).det()).ppartval());
        auto indicator = -heaviside(m_eps - jac(G).det());
        auto E_der = indicator * derJacDet;
//        gsDebugVar(m_evaluator.eval(E_der, pt));

//     现在看来可能是这里有问题？？跟Hugo讨论一下！
        m_assembler.initSystem();
        m_assembler.assemble(E_der);
//        gsDebugVar(m_assembler.matrix());
        result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()), m_assembler.rhs().rows());


//        // CHECK INVERSES OF jac(G)
//        auto JTG = (jac(space) % jac(G));//.tr(); // transpose to switch to rowSpace
//        auto JiG = (jac(space) % jac(G).inv());//.tr(); // transpose to switch to rowSpace
//        auto J_prime = jac(G).det() * JiG;
//        auto Ewinslow = (jac(G).tr() * jac(G)).trace() / jac(G).det();
//        auto Ewinslow_der1 = (2.0 / jac(G).det()) * JTG;//
//        auto Ewinslow_der2 = -(Ewinslow.val() / jac(G).det()) * J_prime;
//        auto Ewinslow_der = Ewinslow_der1 + Ewinslow_der2;
//        auto Euniform_der = 2.0 / pow(S, 2) * (jac(G).det().val() - S.val()) * J_prime;
//
//        // To evaluate expressions at the point 'pt', use the following
////        gsVector<T> pt(2);
////        pt.setConstant(0.5);
////        gsDebugVar(m_evaluator.eval(G, pt));
////        gsDebugVar(m_evaluator.eval(jac(space), pt));
////        gsDebugVar(m_evaluator.eval(jac(G), pt));
////        gsDebugVar(m_evaluator.eval(JTG, pt));
////        gsDebugVar(JTG.Space);
////        gsDebugVar(m_evaluator.eval(JiG, pt));
////        gsDebugVar(JiG.Space);
////
////        gsDebugVar(m_evaluator.eval(J_prime, pt));
////
////        gsDebugVar(m_evaluator.eval(Ewinslow, pt));
////        gsDebugVar(m_evaluator.eval(Ewinslow_der1, pt));
////        gsDebugVar(m_evaluator.eval(Ewinslow_der2, pt));
////        gsDebugVar(m_evaluator.eval(Ewinslow_der, pt));
////        gsDebugVar(m_evaluator.eval(Euniform_der, pt));
//
//        m_assembler.initSystem();
//        m_assembler.assemble(m_lambda1 * Ewinslow_der + m_lambda2 * Euniform_der);
//        result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()), m_assembler.rhs().rows());
    }

protected:
    void convert_design_to_mp(const gsAsConstVector<T> &u, gsMultiPatch<T> &mp) const {
        // Make the geometry
        index_t idx;
        for (size_t i = 0; i != mp.nPatches(); i++) {
            for (size_t c = 0; c != mp.targetDim(); c++) {
                for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)

                    // if it is possible to just loop over the free index
                    // since this function is called very often during optimization
                    if (m_mapper.is_free(k, i, c)) {
                        idx = m_mapper.index(k, i, c);
                        mp.patch(i).coefs()(k, c) = u[idx];
                    }
            }
        }
    }

    gsVector<T> convert_mp_to_gsvec(const gsMultiPatch<T> &mp) const {
        // TODO: Now, it just set all free variables into the vector u
        // I will integrate it with the initialization step

        gsVector<T> currentDesign(m_mapper.freeSize());

        index_t idx;
        for (size_t i = 0; i != mp.nPatches(); i++) {
            for (size_t c = 0; c != mp.targetDim(); c++) {
                for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)
                    // if it is possible to just loop over the free index
                    // since this function is called very often during optimization
                    if (m_mapper.is_free(k, i, c)) {
                        idx = m_mapper.index(k, i, c);
                        currentDesign(idx) = mp.patch(i).coefs()(k, c);
                    }
            }
        }

        return currentDesign;
    }

protected:
    gsMultiPatch<T> m_mp;
    gsDofMapper m_mapper;
    index_t m_patchID;

    gsMultiBasis<T> m_mb;

    mutable gsExprEvaluator<T> m_evaluator;
    mutable gsExprAssembler<T> m_assembler;

    gsOptionList m_options;
    T m_lambda1, m_lambda2;


    /// Number of design variables
    using Base::m_numDesignVars;

    /// Number of constraints
    using Base::m_numConstraints;

    /// Number of nonzero entries in the Constraint Jacobian
    using Base::m_numConJacNonZero;

    /// Lower bounds for the design variables
    using Base::m_desLowerBounds;

    /// Upper bounds for the design variables
    using Base::m_desUpperBounds;

    /// Lower bounds for the constraints
    using Base::m_conLowerBounds;

    /// Upper bounds for the constraints
    using Base::m_conUpperBounds;

    /// Constraint Jacobian non-zero entries rows
    using Base::m_conJacRows;

    /// Constraint Jacobian non-zero entries columns
    using Base::m_conJacCols;

    /// Current design variables (and starting point )
    using Base::m_curDesign;

    mutable gsMatrix<T> m_allBasisVal;
    mutable gsMatrix<T> m_allBasis1stDervs;
    mutable gsMatrix<T> m_allBasis2ndDervs;
    mutable gsMatrix<T> m_gaussWts;
    mutable gsMatrix<index_t> m_gaussIdx;

    T m_eps; // need to handle later, set m_eps = 0.05*S
};


}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBarrierPatch.hpp)
#endif
