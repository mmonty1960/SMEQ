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
#ifndef SMEQ_H
#define SMEQ_H
 
#include "ui_SMEQ.h"
#include <cminpack.h>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_series_data.h>
#include <gsl/gsl_math.h>
#include <qwt_plot_picker.h>
 
class SMEQ : public QWidget, private Ui::SMEQ_DLG
{
    Q_OBJECT
 
public:
    SMEQ(QWidget *parent= nullptr);
    QMap<QString, QLineEdit*> idToLineEdit;
    void setWin(const std::string& _winname);

private:
    Ui::SMEQ_DLG *ui;
    QwtPlot* G1_Rh = nullptr;  // inizializzazione inline, C++11
    QwtPlot* G2_Rns = nullptr;  // inizializzazione inline, C++11
    QwtPlot* Gph = nullptr;  // inizializzazione inline, C++11
    QwtPlotPicker* m_pickerG1 = nullptr;  // cursore con coordinate per G1
    QwtPlotPicker* m_pickerG2 = nullptr;  // cursore con coordinate per G2
    QwtPlotPicker* m_pickerGph = nullptr;  // cursore con coordinate per Gph


 
public Q_SLOTS:
    void setRange();
    void getFileRh();
    void getFileRns();
    void getFileRhSoiled();
    void getFileRnsSoiled();
    double MeanComputing(int iMis);
    void getFileNKover();
    void getFileNKmetal();
    void EMA();
    void EMAcheck();
    void closeEvent ( QCloseEvent * event );
    void InfoFit(int info);
    void PlotRh();
    void PlotRns();
    void PlotRns2();
    void fitRns();
    int nkCheck();
    void calcPlot();
    void bestFit();
    void Gplot(int iGraph, int iRD, int iColor, int iStyle, QString title, int iQ);
    void plotNRloss();
    void RnsPrediction();
    void credits();
}; 
 
#endif
