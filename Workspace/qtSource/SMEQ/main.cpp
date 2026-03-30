/*
Solar Mirror Exhaustive Qualification (SMEQ) software

SMEQ is based on two milestones:
    1) the Equivalent Model Algorithm for Solar Mirror (EMA4SM)
    2) the Soiled Mirror Modeling software (SoilMirMod)

EMA4SM allows to model the hemispherical spectral reflectance
    and predicts the near-specular reflectance at any incidence angle
    and any acceptance angle in the range covered by suitable experimental data
    such as SMQ2, VLABS, UNIZA and S2R at some spare wavelengths
    [AIP Conference Proceedings 2303, 100005 (2020)]
    DOI: https://doi.org/10.1063/5.0028768

SoilMirMod software deals with the modeling of the reflectance spectral loss
    due to soiling. The reflectance loss is modeled by 3 different phenomema
    induced by the soiling particles:
    i) absorption; ii) scattering; iii) diffraction
    [Solar Energy 259 (2023) 356-363]
    DOI: https://doi.org/10.1016/j.solener.2023.05.017

Author: Marco Montecchi
        ENEA-Casaccia
        marco.montecchi@enea.it

This file is part of SMEQ.

    SMEQ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation version 3 of the License

    SMEQ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2026 Marco Montecchi
*/
#include <QApplication>
#include "SMEQ.h"
 
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    SMEQ *dialog = new SMEQ;
    dialog->show();
    return app.exec();
}

