/** @file gsBarrierPatch.hpp

    @brief Provides Coons's patch construction from boundary data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Haberleitner, A. Mantzaflaris
*/

#pragma once

#include<gsModeling/gsBarrierPatch.h>
//#include<gsNurbs/gsTensorBSplineBasis.h>
//#include<gsNurbs/gsBSpline.h>
//#include<gsTensor/gsGridIterator.h>

namespace gismo
{
template<short_t d, typename T>
gsBarrierPatch<d,T>::gsBarrierPatch(const gsMultiPatch<T> &bRep,
                const index_t &method, const bool &plot_init)
        :
//m_mapper(mapper),
        m_bRep(bRep),
        //m_mp(mp),
        m_method(method),
        m_plot(false),
        m_plot_init(plot_init)
{
    GISMO_ASSERT(d == 2, "This method is only available for d==2");
    // TODO: assign m_mp by using initialization
    // m_mp = ;

    //sanityCheckInput

//        m_eps = 0.0 * m_area;

    ///*gsDebug << &m_bRep.basis(0)->knots() << "\n";*/

    //gsBasis<T>* westBoundary = &m_bRep.basis(0);
    //gsDebug << typeid(westBoundary).name() << "\n";

    //gsDebug << westBoundary->knots() << "\n";

    /*gsNurbsBasis<T> westBoundary = static_cast<gsNurbsBasis<T>&> (m_bRep.basis(0));
    gsDebug << typeid(westBoundary).name() << "\n";
    gsDebug << westBoundary.knots() << "\n";*/

    initialization();
    m_mb = gsMultiBasis<T>(m_mp);
    makeMapper();

    m_area = computeArea();
    m_eps = 0.05 * m_area; // ATTENTION!!

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
    m_curDesign = convert_mp_to_gsvec();
}

template <short_t d, typename T>
const gsGeometry<T> & gsBarrierPatch<d,T>::compute()
{
    GISMO_NO_IMPLEMENTATION;
}

template <short_t d, typename T>
gsObjQualityImprovePt<d,T>::gsObjQualityImprovePt(const gsMultiPatch<T> &patches,
                                                  const gsDofMapper &mapper,
                                                  const T eps)
:
m_mp(patches),
m_mapper(mapper),
m_mb(m_mp),
m_patchID(0),
m_eps(eps)
{

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

template <short_t d, typename T>
T gsObjQualityImprovePt<d,T>::evalObj(const gsAsConstVector<T> &u) const
{
    T F = 0;
    gsMultiPatch<T> mp = m_mp;
    this->convert_design_to_mp(u, mp);

    geometryMap G = m_evaluator.getMap(mp);

    //ReLU, Rectified Linear Unit, how to implement?

//        auto jacGdet = jac(G).det();
//        auto EfoldElimination = jacGdet;//.val();
//        auto EfoldElimination = (m_eps-jac(G).det()).ppartval();//.val();

    gsConstantFunction<T> lambda1(m_lambda1, d);
    gsConstantFunction<T> lambda2(m_lambda2, d);
    auto lam1 = m_evaluator.getVariable(lambda1);
    auto lam2 = m_evaluator.getVariable(lambda2);

    gsConstantFunction<T> penaltyFactor(1e20, d);
    auto penaF = m_evaluator.getVariable(penaltyFactor);

    auto Ewinslow = (jac(G) % jac(G)) / jac(G).det();
//        auto Ewinslow2 = (jac(G).tr() * jac(G)).trace() / jac(G).det();
    auto Euniform = pow(jac(G).det(), 2);
//        gsDebugVar(m_eps);
    auto EimprovePtReal = lam1 * Ewinslow + lam2 * Euniform;


//        if jac(G).det() < 0, then indicator = 1; otherwise indicator = 0
//        auto indicator = (-jac(G).det()).ppartval()/(1e-20-(-(jac(G).det()).ppartval()));
//        auto indicator2 = -jac(G).det().ppartval()/(1e-20-((jac(G).det()).ppartval()));

    auto indicator = heaviside(jac(G).det());
    auto indicatorNot = heaviside(-jac(G).det());
    auto EimprovePt = indicator * EimprovePtReal + indicatorNot * penaF;

    // 这里需要一个指示函数，目前只对不为0的部分做测试
//        gsConstantFunction<T> lambda1(-0.5, d);
//        auto lam1 = m_evaluator.getVariable(lambda1);
//        auto EfoldElimination = jac(G).det();
//        gsVector<T> pt(2);
//        pt.setConstant(0.3);
//        gsDebugVar(m_evaluator.eval(indicator, pt));
//        gsDebugVar(m_evaluator.eval(indicator2, pt));
//        gsDebugVar(m_evaluator.eval(EimprovePt, pt));
//
//        gsDebugVar(m_evaluator.eval(Ewinslow2, pt));
//        gsDebugVar(m_evaluator.eval(m_eps-jac(G).det(), pt));

//        auto EfoldElimination = jac(G).det() * jac(G).det(); // for test purpose
    F = m_evaluator.integral(EimprovePt);

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

}// namespace gismo
