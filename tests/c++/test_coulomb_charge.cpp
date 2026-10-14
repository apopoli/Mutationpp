/*
 * Copyright 2026 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "mutation++.h"
#include "CoulombIntegrals.h"
#include <catch.hpp>

using namespace Mutation::Transport;
using namespace Catch;

/**
 * For the screened Coulomb potential z1 z2 e^2/(4 pi eps0 r) exp(-r/lambda),
 * b = |z1 z2| e^2/(8 pi eps0 k T) depends on T and z1 z2 only through
 * T/|z1 z2|. At fixed Debye length the integrals for a charge product zz at
 * temperature T must therefore equal the singly charged ones at T/zz.
 */
TEST_CASE("Debye-Huckel integrals scale with the charge product",
    "[transport]")
{
    const CoulombType types[] = {
        Q11_ATT, Q11_REP, Q22_ATT, Q22_REP, Q14_ATT, Q15_REP, Q24_ATT,
        BST_ATT, CST_REP, EST_ATT, Q12_ATT, Q13_REP, Q23_ATT, AST_REP
    };
    const double charge_products[] = {2.0, 3.0, 4.0, 6.0, 9.0};
    const double temperatures[] = {3000.0, 10000.0, 20000.0, 40000.0};

    for (double zz : charge_products) {
        DebyeHuckleEvaluator single, multiple;
        multiple.setChargeProduct(zz);
        CHECK(single.chargeProduct() == 1.0);
        CHECK(multiple.chargeProduct() == zz);

        // Same Debye length for both evaluators
        single.setDebyeLength(10000.0, 1.0e22);
        multiple.setDebyeLength(10000.0, 1.0e22);

        for (double T : temperatures) {
            for (CoulombType type : types) {
                const double ref = single(T / zz, type);
                CHECK(multiple(T, type) == Approx(ref).epsilon(1.0e-12));
            }
        }
    }
}

/**
 * The integrals grow with the charge product: at fixed T and lambda a doubly
 * charged pair scatters more than a singly charged one.
 */
TEST_CASE("Debye-Huckel integrals increase with the charge product",
    "[transport]")
{
    DebyeHuckleEvaluator e1, e2, e4;
    e2.setChargeProduct(2.0);
    e4.setChargeProduct(4.0);
    for (auto* e : {&e1, &e2, &e4})
        e->setDebyeLength(20000.0, 1.0e23);

    for (CoulombType type : {Q11_ATT, Q11_REP, Q22_ATT, Q22_REP}) {
        const double q1 = e1(20000.0, type);
        const double q2 = e2(20000.0, type);
        const double q4 = e4(20000.0, type);
        CHECK(q2 > q1);
        CHECK(q4 > q2);
    }
}
