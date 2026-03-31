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
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>
#include <qtextstream.h>
#include "SMEQ.h"
#include <cminpack.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_series_data.h>
#include <qwt_plot_grid.h>
#include <qwt_legend.h>
#include <QwtPickerMachine>
#include <QwtPlotZoomer>
#include <qwt_text.h>
#include <qwt_symbol.h>
#include <gsl/gsl_math.h>
#include <tgmath.h>
#include <QDesktopServices>

using namespace std;

//global variables *******************************
QColor myColor[7]={Qt::black,Qt::blue,Qt::cyan,Qt::green,Qt::magenta,Qt::red,Qt::yellow};
QString pathroot;
int Nw=437;
int iLock=0;
static double xy[2][100];
static double SUMw,d_layer,Nini,Kini,thetaService,phiService,TISc,Pa[5];
static double Pig=acos(-1.);
static double REMAs,REMAp,REMAu,RnsEMAs,RnsEMAp,RnsEMAu,R1s,R1p,nOl,kOl;
double ParFit[8][2];
static double MRh[1000][30];
/*
  MRh[1000][30] main data storage matrix
        [i][0] = wavelength (nm)
        [i][1] = weight (absolute)
        [i][2] = R_reference
        [i][3] = n over-layer
        [i][4] = k over-layer
        [i][5] = n metallic-layer
        [i][6] = k metallic-layer
        [i][7] = Rs simulation
        [i][8] = Rp simulation
        [i][9] = Runpolarised simulation
        [i][10] = Rh_clean experimental
        [i][11] = Rh_clean simulation
        [i][12] = Rns_clean experimental
        [i][13] = Rns_clean simulation
        ...
        [i][20] = Rh_soil experimental
        [i][21] = Rh_soil simulation
        [i][22] = Rns_soil experimental
        [i][23] = Rns_soil simulation
        [i][24]= (RhSoil-RhClean)/RhClean exp
        [i][25]= (RhSoil-RhClean)/RhClean calc
        [i][26]= (RnsSoil-RnsClean)/RnsClean exp
        [i][27]= (RnsSoil-RnsClean)/RnsClean calc
        [i][28]= RnsSoiledPrediction
*/

// funzioni invocate
void parabolicFit(int n,double& A, double& B, double& C);
void spada(int column, QString fileRh);
void RsRpEMA(int iWL, double theta, double phiS, int iTISopt);
void RsRp(complex<double> N1, complex<double> N2, complex<double> theta1,
          double& Rs,double& Rp, complex<double>& theta2);
void FindRoot(int ix,int i,int iWarming);
double FT(int ix,int i,double t);
double FunFit(double xp);
double TIScomplete(double phiS, double phiP,double lambdaOL, double thetaOL);
double fmodel(double wl,double theta,double A,double S,double sigma,double C,double D,double k,double L,double B,double E);
int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag);//fit SoilMirModel
int fcn0(void *p, int m, int n, const double *x, double *fvec, int iflag);//fit Rnearspecular with EMA4SM

struct pointToFit2{
    int NdatFit;
    int Nmin;
    int Nmax;
}pTF2[1];



SMEQ::SMEQ(QWidget *parent){
    setupUi(this); // this sets up GUI
    Q_UNUSED( parent )

    // signals/slots mechanism in action
    connect(comboBox_standard,    SIGNAL(currentIndexChanged(int)), this, SLOT(setRange()));
    connect(pushButton_loadRh,    SIGNAL(clicked()), this, SLOT( getFileRh() )     );
    connect(pushButton_loadRns,   SIGNAL(clicked()), this, SLOT( getFileRns() )    );
    connect(pushButton_loadRhSoiled,    SIGNAL(clicked()), this, SLOT( getFileRhSoiled() )     );
    connect(pushButton_loadRnsSoiled,   SIGNAL(clicked()), this, SLOT( getFileRnsSoiled() )    );
    connect(pushButton_PlotRns,   SIGNAL(clicked()), this, SLOT( PlotRns() )       );
    connect(comboBox_nkOL,        SIGNAL(currentIndexChanged(int)), this, SLOT( getFileNKover() ) );
    connect(comboBox_nkMetal,     SIGNAL(currentIndexChanged(int)), this, SLOT( getFileNKmetal()) );
    connect(pushButton_EMA,       SIGNAL(clicked()), this, SLOT( EMA() )           );
    connect(pushButton_check,     SIGNAL(clicked()), this, SLOT( EMAcheck() )      );
    connect(pushButton_fitRns,    SIGNAL(clicked()), this, SLOT( fitRns() )        );
    connect(pushButton_plotNormLoss,SIGNAL(clicked()), this, SLOT(plotNRloss() )   );
    connect(pushButton_BF,SIGNAL(clicked()), this, SLOT(bestFit() )   );
    connect(pushButton_RnsPrediction,SIGNAL(clicked()), this, SLOT(RnsPrediction() )   );
    connect(pushButton_plotRh,    SIGNAL(clicked()), this, SLOT(PlotRh() )   );
    connect(pushButton_PlotRns2,  SIGNAL(clicked()), this, SLOT(PlotRns2() )   );
    connect(pushButton_credits,   SIGNAL(clicked()), this, SLOT(credits() )   );

#ifdef __unix__
#define IS_POSIX 1
#else
#define IS_POSIX 0
#endif

    QDir dir;  //current directory
    dir.cdUp();//cd ..
    dir.cdUp();//cd ..
    if (IS_POSIX == 1) {
        //Linux path initialization
        //nothing to do
    }
    else {
        //windows path inizialization
        dir.cdUp();//cd ..
    }
    pathroot=dir.absolutePath()+"/SMEQ";

    setRange();//initialization

//    idToCheckBox["checkBox_Pa0"]=checkBox_Pa0;
//    idToCheckBox["checkBox_Pa1"]=checkBox_Pa1;
//    idToCheckBox["checkBox_Pa2"]=checkBox_Pa2;
//    idToCheckBox["checkBox_Pa3"]=checkBox_Pa3;
//    idToCheckBox["checkBox_Pa4"]=checkBox_Pa4;
}

void SMEQ::closeEvent ( QCloseEvent *  ){
    qApp->quit();
}

void SMEQ::credits(){
    QString pdfPath = QDir(QCoreApplication::applicationDirPath())
    .absoluteFilePath("SMEQ_Credits.pdf");
    QDesktopServices::openUrl(QUrl::fromLocalFile(pdfPath));
}


void SMEQ::setRange(){
    int iStd = comboBox_standard -> currentIndex();
    comboBox_range -> clear();
    comboBox_range -> addItem("solar 320-2500 nm");
    comboBox_range -> addItem("UVA 315-400 nm");
    if(iStd==0 || iStd==1){
        comboBox_range -> addItem("UVB 280-315 nm");
        comboBox_range -> addItem("UV  280-400 nm");
    }
}


double SMEQ::MeanComputing(int iMis){
    //computing mean value with the selected standard spectrum
    double SUMwR=0.;
    for(int i=0;i<Nw;i++){
        SUMwR=SUMwR+MRh[i][1]*MRh[i][iMis];
    }
    double meanVal=SUMwR/SUMw;
    printf("->MeanComputig: iMis=%d Ndat=%d SUMw=%f SUMwR=%f\n",iMis,Nw,SUMw,SUMwR);
    return(meanVal);
}


void SMEQ::Gplot(int iGraph, int iRD, int iColor, int iStyle, QString title, int iQ){
    //int iGraph= 1: G1_Rh
    //            2: G2_Rns
    //            3: Gph
    //int iRD= 0->redraw
    //         1->overwrite
    //int iColor= color index 0->black
    //                        1->blue
    //                        2->cyan
    //                        3->green
    //                        4->magenta
    //                        5->red
    //                        6->yellow
    //int iStyle: 0->SolidLine
    //            1->DotLine
    //            2->DashLine

    printf("->Gplot: iGraph=%d iColor=%d iStyle=%d iQuantity=%d Ndata=%d\n",iGraph,iColor,iStyle,iQ,Nw);

    vector <double> Xp(Nw),Yp(Nw);
    for(int ii=0;ii< Nw;ii++){
        Xp[ii]=MRh[ii][0]; //WL
        Yp[ii]=MRh[ii][iQ];//quantity to be plot
    }
    // Make the grid on
    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->setPen(QPen(Qt::gray, 0.0, Qt::DotLine));
    grid->enableX(true);
    grid->enableXMin(true);
    grid->enableY(true);
    grid->enableYMin(true);
    QwtPlotCurve *curve1=new QwtPlotCurve("Curve 1");
    curve1->setSamples(Xp.data(), Yp.data(), Nw);
    if(iStyle==0)
        curve1->setPen(QPen(myColor[iColor],3,Qt::SolidLine));
    else if(iStyle==1)
        curve1->setPen(QPen(myColor[iColor],3,Qt::DotLine));
    else
        curve1->setPen(QPen(myColor[iColor],3,Qt::DashLine));
    if(iGraph==1){
        if(G1_Rh == nullptr) {
            //the window does not exist
            G1_Rh = new QwtPlot();
            connect(G1_Rh, &QObject::destroyed, this, [this]() {
                G1_Rh = nullptr;
                m_pickerG1 = nullptr;
            });
            G1_Rh->setAttribute(Qt::WA_DeleteOnClose);
            m_pickerG1 = new QwtPlotPicker(
                QwtAxis::XBottom, QwtAxis::YLeft,
                QwtPlotPicker::CrossRubberBand,
                QwtPicker::AlwaysOn,
                G1_Rh->canvas()
                );
            m_pickerG1->setStateMachine(new QwtPickerDragPointMachine());
            m_pickerG1->setRubberBandPen(QPen(Qt::red));
            m_pickerG1->setTrackerPen(QPen(Qt::black));
        }
        else {
            // esiste già: portalo in primo piano
            G1_Rh->raise();
            G1_Rh->activateWindow();
        }
        if(iRD==0)
            G1_Rh->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
        G1_Rh -> setTitle(title);
        G1_Rh -> setAxisTitle(0,"Rhemispherical");
        G1_Rh -> setAxisTitle(2,"wavelength (nm)");
        G1_Rh -> setAxisScale(2,Xp[0],Xp[Nw-1],0.);
        curve1->attach(G1_Rh);
        grid->attach(G1_Rh);
        G1_Rh->replot();
        G1_Rh->show();
    }
    else if(iGraph==2){
        if(G2_Rns == nullptr) {
            //the window does not exist
            G2_Rns = new QwtPlot();
            connect(G2_Rns, &QObject::destroyed, this, [this]() {
                G2_Rns = nullptr;
                m_pickerG2 = nullptr;
            });
            G2_Rns->setAttribute(Qt::WA_DeleteOnClose);
            m_pickerG2 = new QwtPlotPicker(
                QwtAxis::XBottom, QwtAxis::YLeft,
                QwtPlotPicker::CrossRubberBand,
                QwtPicker::AlwaysOn,
                G2_Rns->canvas()
                );
            m_pickerG2->setStateMachine(new QwtPickerDragPointMachine());
            m_pickerG2->setRubberBandPen(QPen(Qt::red));
            m_pickerG2->setTrackerPen(QPen(Qt::black));
        }
        else {
            // esiste già: portalo in primo piano
            G2_Rns->raise();
            G2_Rns->activateWindow();
        }
        if(iRD==0)
            G2_Rns->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
        G2_Rns -> setTitle("R-nearspecular spectrum");
        G2_Rns -> setAxisTitle(0,"Reflectance");
        G2_Rns -> setAxisTitle(2,"Wavelength (nm)");
        G2_Rns -> setAxisScale(2,Xp[0],Xp[Nw-1],0.);
        curve1->attach(G2_Rns);
        grid->attach(G2_Rns);
        G2_Rns->replot();
        G2_Rns->show();
    }
    else if(iGraph==3){
        if (Gph == nullptr) {
            //the window does not exist
            Gph = new QwtPlot();
            connect(Gph, &QObject::destroyed, this, [this]() {
                Gph = nullptr;
                m_pickerGph = nullptr;
            });
            Gph->setAttribute(Qt::WA_DeleteOnClose);
            m_pickerGph = new QwtPlotPicker(
                QwtAxis::XBottom, QwtAxis::YLeft,
                QwtPlotPicker::CrossRubberBand,
                QwtPicker::AlwaysOn,
                Gph->canvas()
                );
            m_pickerGph->setStateMachine(new QwtPickerDragPointMachine());
            m_pickerGph->setRubberBandPen(QPen(Qt::red));
            m_pickerGph->setTrackerPen(QPen(Qt::black));
        }
        else {
            // esiste già: portalo in primo piano
            Gph->raise();
            Gph->activateWindow();
        }
        if(iRD==0)
            Gph->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
        double Ymax=doubleSpinBox_Ymax->value();
        double Ymin=doubleSpinBox_Ymin->value();
        Gph -> setTitle("Normalized-reflectance loss");
        Gph -> setAxisTitle(0,"norm. loss");
        Gph -> setAxisTitle(2,"Wavelength (nm)");
        Gph -> setAxisScale(2,Xp[0],Xp[Nw-1],0.);
        Gph -> setAxisScale(0,Ymin,Ymax);
        curve1->attach(Gph);
        grid->attach(Gph);
        Gph->replot();
        Gph->show();
    }
}

void SMEQ::getFileRh(){
    int iRg,iStd,iRowB=0,iRowE=0;
    QString fileRh,fileW;

    iRg=comboBox_range -> currentIndex();     //range
    iStd=comboBox_standard -> currentIndex(); //standard
    if(iStd==0)
        fileW=pathroot+"/Standard/IEC60904b3Step5nm.txt";
    else if(iStd==1)
        fileW=pathroot+"/Standard/ASTMG173SP.txt";
    else if(iStd==2)
        fileW=pathroot+"/Standard/ISO9050.txt";
    else if(iStd==3)
        fileW=pathroot+"/Standard/ISO9845b1.txt";
    else if(iStd==4)
        fileW=pathroot+"/Standard/E891.txt";
    if(iStd<=1){
        if(iRg==0){
            iRowB=9;
            iRowE=445;
        }
        else if(iRg==1){
            iRowB=8;
            iRowE=25;
        }
        else if(iRg==2){
            iRowB=1;
            iRowE=8;
        }
        else if(iRg==3){
            iRowB=1;
            iRowE=25;
        }
    }
    else{
        if(iRg==0){
            iRowB=5;
            iRowE=441;
        }
        else if(iRg==1){
            iRowB=4;
            iRowE=21;
        }
    }
    QFile fW(fileW);
    if (!fW.open(QIODevice::ReadOnly | QIODevice::Text)){
        printf("ERROR opening file= %s\n",fileW.toStdString().c_str());
        return;
    }
    QTextStream stream0 (&fW);
    Nw=0;
    SUMw=0.;
    int iRow=0;
    double dum1,dum2;
    stream0.readLine();
    stream0.readLine();
    do {
        iRow++;
        if(iRowB<=iRow && iRow<=iRowE){
            stream0 >> MRh[Nw][0] >> MRh[Nw][1];
            MRh[Nw][0]=MRh[Nw][0]/10.;//Angtstrom -> nm
            //printf("%d	%f	%f \n", Nw,MRh[Nw][0],MRh[Nw][1]);
            SUMw=SUMw+MRh[Nw][1];
            Nw++;
        }
        else {
            stream0 >> dum1 >> dum2;
        }
    } while (!stream0.atEnd());
    printf("Loaded %s \n Nw=%d SUMw=%f\n",(fileW.toStdString()).c_str(),Nw,SUMw);
    fW.close();
    lineEdit_Nwl->setText("Nwl= "+QString::number(Nw));

    //ex R_reference
    for(int i=0;i<Nw;i++) MRh[i][2]=1.0;
    
    //Rh file
    fileRh = QFileDialog::getOpenFileName(
                this,
                "Choose a Rhemispherical exp. file",
                pathroot+"/Spectra");
    if(fileRh.isEmpty())
        return;
    lineEditRh -> setText(fileRh.section('/',-1));
    spada(10,fileRh);
    printf("Loaded %s\n",(fileRh.toStdString()).c_str());
    double meanVal=MeanComputing(10);
    lineEdit_meanRh->setText(QString::number(meanVal,'f',4));
    PlotRh();
}

void SMEQ::PlotRh(){
    //Plot Rh exp
    Gplot(1,0,0,0,"R-hemispherical spectrum",10);
}


int SMEQ::nkCheck(){
    int iNKol=comboBox_nkOL->currentIndex();
    int iNKmetal=comboBox_nkMetal->currentIndex();
    if(iNKol<0 || iNKmetal<0){
        QMessageBox msgBox;
        msgBox.setText("Please select nk ov over-layer and metallic-layer");
        msgBox.exec();
        return(-1);
    }
    return(0);
}


void SMEQ::getFileNKover(){
    QString fileNK;
    int iOL=comboBox_nkOL->currentIndex();
    if(iOL==0)
        fileNK=pathroot+"/FileNK/GlassBK7_nk.txt";
    else
        fileNK=pathroot+"/FileNK/sio2_nk.txt";
    spada(3,fileNK);
    printf("Loaded %s\n",(fileNK.toStdString()).c_str());
}



void SMEQ::getFileNKmetal(){
    QString fileNK;
    int iMetal=comboBox_nkMetal->currentIndex();
    if(iMetal==0)
        fileNK=pathroot+"/FileNK/Silver_nk.txt";
    else
        fileNK=pathroot+"/FileNK/Aluminum_nk.txt";
    spada(5,fileNK);
    printf("Loaded %s\n",(fileNK.toStdString()).c_str());
}



void SMEQ::EMAcheck(){
    int iOK=nkCheck();
    if(iOK<0)
        return;
    double term,Drms,Dmin,Dmax;
    double theta=dSB_thetaRh->value();
    theta=theta/180.*Pig;

    d_layer=dSBdoverlayer -> value();
    d_layer=d_layer*1.e+06;//mm -> nm

    Drms=0.;
    Dmin=1.e+36;
    Dmax=-1.e+36;

    //PLOT Rh experimental
    Gplot(1,0,0,0,"check EMA",10);

    //Plot Rh EMA
    for(int ii=0;ii< Nw;ii++){
        RsRpEMA(ii,theta,0.,-1);
        MRh[ii][11]=REMAu;
        term=REMAu-MRh[ii][10];
        Dmin=min(term,Dmin);
        Dmax=max(term,Dmax);
        Drms=Drms+term*term;
    }
    Gplot(1,1,3,2,"check EMA",11);
    lineEdit_RMS -> setText(QString::number(sqrt(Drms/double(Nw))));
    lineEdit_PV -> setText(QString::number(Dmax-Dmin));
    double meanVal=MeanComputing(11);
    lineEdit_meanRhEMA->setText(QString::number(meanVal,'f',4));
}




void SMEQ::EMA(){
    int iOK=nkCheck();
    if(iOK<0)
        return;
    int ix,iAlarm=0;
    double SUMwRh,term,Drms,Dmin,Dmax;
    double thetaRh=dSB_thetaRh ->value();
    thetaRh=thetaRh/180.*Pig;
    thetaService=thetaRh;//to be pass to FT
    QString fileMatNK=pathroot+"/MatNK.txt";
    QFile fileMnk(fileMatNK);
    fileMnk.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream streamMnk ( &fileMnk );

    d_layer=dSBdoverlayer -> value();
    d_layer=d_layer*1.e+06;//mm -> nm

    ix=comboBox_ix -> currentIndex();
    if(ix==0)
        ix=4;//k over
    else if(ix==1)
        ix=5;//n metal
    else if(ix==2)
        ix=6;//k metal
    else if(ix==3)
        ix=56;//n&k metal
    Drms=0.;
    Dmin=1.e+36;
    Dmax=-1.e+36;

    //PLOT Rh experimental
    Gplot(1,0,0,0,"EMA computing",10);

    // Plot initial simulated Rh
    for(int ii=0;ii< Nw;ii++){
        RsRpEMA(ii,thetaRh,0.,-1);
        MRh[ii][11]=REMAu;
    }
    Gplot(1,1,3,2,"EMA computing",11);

    SUMwRh=0.;
    for(int i=0;i<Nw;i++){
        FindRoot(ix,i,0);
        RsRpEMA(i,thetaRh,0.,-1);
        //    printf("i=%d WL=%f Rs=%f Rp=%f\n",i,MRh[i][0],REMAs,REMAp);
        MRh[i][7]=REMAs;
        MRh[i][8]=REMAp;
        MRh[i][9]=REMAu;
        term=MRh[i][10]*MRh[i][1];
        Dmin=min(REMAu-MRh[i][10],Dmin);
        Dmax=max(REMAu-MRh[i][10],Dmax);
        Drms=Drms+(REMAu-MRh[i][10])*(REMAu-MRh[i][10]);
        SUMwRh=SUMwRh+term;
        if(iAlarm==0 && REMAs*REMAp<0.){
            QMessageBox msgBox;
            msgBox.setText("Rs or/and Rp <0");
            msgBox.exec();
            iAlarm=1;
        }
        streamMnk << MRh[i][0] <<"\t"<< MRh[i][3] <<"\t"<< MRh[i][4] <<"\t"<< MRh[i][5] <<"\t"<< MRh[i][6] <<"\n";
    }
    fileMnk.close();
    lineEdit_RMS -> setText(QString::number(sqrt(Drms/double(Nw))));
    lineEdit_PV -> setText(QString::number(Dmax-Dmin));
    //PLOT final simulated Rh
    for(int ii=0;ii< Nw;ii++){
        RsRpEMA(ii,thetaRh,0.,-1);
        MRh[ii][11]=REMAu;
    }
    Gplot(1,1,5,2,"EMA computing",11);
}


void SMEQ::getFileRns(){
    QString fileRns = QFileDialog::getOpenFileName(
        this,
        "Choose a Rnear_specular exp. file",
        pathroot+"/Spectra");
    if(fileRns.isEmpty())
        return;
    lineEdit_Rns -> setText(fileRns.section('/',-1));
    spada(12,fileRns);
    printf("Loaded %s\n",(fileRns.toStdString()).c_str());
    double meanVal=MeanComputing(12);
    lineEdit_meanRnsExp->setText(QString::number(meanVal,'f',4));
    PlotRns();
}


void SMEQ::PlotRns(){
    //plot of Rns experimental and computed
    //plot Rns exp
    Gplot(2,0,0,0,"R-nearspecular spectrum",12);

    //plot Rns computed by EMA4SM
    double theta=dSB_thetaRns->value();
    theta=theta/180.*Pig;
    double phiS=dSB_Phi->value();
    phiS=phiS/1000;
    double sigma=dSB_sigmaClean->value();
    TISc=doubleSpinBox_TIS->value();
    printf("->PlotRns computing with theta(rad)=%f phiS(rad)=%f sigma(nm)=%f\n",theta,phiS,sigma);
    Pa[0]=sigma;//sigma (nm)
    double sum=0.;
    double rms=0;
    for(int ii=0;ii< Nw;ii++){
        RsRpEMA(ii, theta, phiS, 0);
        rms=rms+pow(RnsEMAu-MRh[ii][12],2.);
        MRh[ii][13]=RnsEMAu;
        sum=sum+MRh[ii][1]*MRh[ii][13];
    }
    Gplot(2,1,3,2,"R-nearspecular spectrum",13);
    lineEdit_RMSfit->setText(QString::number(sqrt(rms)/double(Nw),'e',4));
    double meanVal=MeanComputing(13);
    lineEdit_meanRnsCalc->setText(QString::number(meanVal,'f',4));
}



void SMEQ::fitRns(){
    int info;
    int n=1;
    int m=Nw;
    int lwa=m*n+5*n+m;
    vector <int> iwa(n);
    vector <double> x(n),fvec(m),wa(lwa);
    double tol=sqrt(dpmpar(1));
    thetaService=dSB_thetaRns->value();
    thetaService=thetaService/180.*Pig;//to be pass to fcn0
    phiService=dSB_Phi->value();
    phiService=phiService/1000;//to be pass to fcn0
    TISc=doubleSpinBox_TIS->value();
    x[0]=dSB_sigmaClean->value();//sigma
    if(m>0 && n>0){
        printf("Launch bestFit with Np=%d Ndata=%d data theta(rad)=%f phi(rad)=%f sigma(nm)=%f\n",n,m,thetaService,phiService,x[0]);
        info=lmdif1(fcn0, &pTF2, m, n, x.data(),fvec.data(), tol, iwa.data(), wa.data(), lwa);
    }
    else{
        printf("Nmis=0 || m=0 || n=0 then RETURN!");
        return;
    }
    dSB_sigmaClean->setValue(x[0]);
    printf(" done with info=%d\n",info);
    InfoFit(info);
    PlotRns();
}



void SMEQ::InfoFit(int info){
    if(info==0)
        lineEdit_infoFit -> setText("Compute: 0= improper input parameters");
    else if(info==1)
        lineEdit_infoFit -> setText("Compute: 1= relative error in the sum of squares is at most tol");
    else if(info==2)
        lineEdit_infoFit -> setText("Compute: 2= relative error between x and the solution is at most tol");
    else if(info==3)
        lineEdit_infoFit -> setText("Compute: 3= conditions for info = 1 and info = 2 both hold");
    else if(info==4)
        lineEdit_infoFit -> setText("Compute: 4= fvec is orthogonal to the columns of the jacobian to machine precision");
    else if(info==5)
        lineEdit_infoFit -> setText("Compute: 5= number of calls to fcn has reached or exceeded 200*(n+1)");
    else if(info==6)
        lineEdit_infoFit -> setText("Compute: 6= tol is too small. No further reduction in SUM_sqare is possible");
    else if(info==7)
        lineEdit_infoFit -> setText("Compute: 7= tol is too small. no further improvement in the approximate solution x is possible");
}



int fcn0(void *p, int m, int n, const double *x, double *fvec, int iflag){
    //fit Rnearspecular with EMA4SM
    int iv=0;
    Pa[0]=x[0];
    double chi2=0.;
    for(int j=0;j<Nw;j++){
        RsRpEMA(j, thetaService, phiService, 0);
        fvec[iv]=RnsEMAu-MRh[j][12];
        chi2=chi2+fvec[iv]*fvec[iv];
        //printf("RnsEMAu=%f RnsExp=%f fvec[iv]=%f\n",RnsEMAu,MRh[j][12],fvec[iv]);
        iv++;
    }
    printf("\t->fcn0: sigma=%f chi2=%f\n",Pa[0],chi2/Nw);
    return(0);
}

void parabolicFit(int n, double& A, double& B, double& C){
    long double matrix[3][4], ratio, a;
    long double sum_x=0.,sum_y=0.,sum_x2=0.,sum_x3=0.,sum_x4=0.,sum_xy=0.,sum_x2y=0.;
    int i, j , k;
    for(i = 0; i < n; i++){
        sum_x += xy[0][i];
        sum_y += xy[1][i];
        sum_x2 += pow(xy[0][i], 2);
        sum_x3 += pow(xy[0][i], 3);
        sum_x4 += pow(xy[0][i], 4);
        sum_xy += xy[0][i]*xy[1][i];
        sum_x2y += pow(xy[0][i], 2) * xy[1][i];
    }
    matrix[0][0] = n;
    matrix[0][1] = sum_x;
    matrix[0][2] = sum_x2;
    matrix[0][3] = sum_y;
    matrix[1][0] = sum_x;
    matrix[1][1] = sum_x2;
    matrix[1][2] = sum_x3;
    matrix[1][3] = sum_xy;
    matrix[2][0] = sum_x2;
    matrix[2][1] = sum_x3;
    matrix[2][2] = sum_x4;
    matrix[2][3] = sum_x2y;
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            if(i!=j){
                ratio = matrix[j][i]/matrix[i][i];
                for(k = 0; k < 4; k++){
                    matrix[j][k] -= ratio * matrix[i][k];
                }
            }
        }
    }
    for(i = 0; i < 3; i++){
        a = matrix[i][i];
        for(j = 0; j < 4; j++){
            matrix[i][j] /= a;
        }
    }
    C=matrix[0][3];
    B=matrix[1][3];
    A=matrix[2][3];
}


void spada(int ic, QString fileRh){
    int i,Ndata,Nflag,Nflag0,n,istep,iAlarmMIN=0,iAlarmMAX=0;
    double A,B,C;
    static double RhExp[5000][5];
    double Wa,Wb,WL,Yfin,DET;//err
    QString line,msg;
    QStringList list;
    printf("spada ic=%d file=%s",ic,fileRh.toStdString().c_str());
    
    QFile fRh(fileRh);
    if (!fRh.open(QIODevice::ReadOnly | QIODevice::Text)){
        printf("ERROR opening %s\n",fileRh.toStdString().c_str());
        return;
    }
    QTextStream stream ( &fRh );
    Ndata=0;
    double val[5];
    int Narg=2;
    if(ic == 3 || ic == 5)
        Narg=5;
    int Nskip=0;
    printf("\tNarg=%d\n",Narg);
    do{
        line = fRh.readLine();
        line=line.simplified();
        list = line.split(QRegularExpression("\\s+"));
        int nV=list.size();
        bool OK=false;
        if(nV==Narg){
            for(int index=0;index<nV;index++){
                val[index]=list.at(index).toDouble(&OK);
                if(!OK)
                    index=nV;
            }
            if(OK){
                for (int k=0;k<nV;k++){
                    RhExp[Ndata][k]=val[k];
                    //cout<<val[k]<<"\t";
                }
                //cout<<"\n";
                Ndata++;
            }
            else
                Nskip++;
        }
        else{
            Nskip++;
            //cout<<"   ->>>>> "<<nV<<"\t"<<OK<<"\n";
        }
    }while(!stream.atEnd());

    printf("Ndati=%d Nskip_line=%d\n",Ndata,Nskip);
    fRh.close();
    if(RhExp[0][0] < RhExp[Ndata-1][0]){
        istep=1;
        Nflag0=0;
        if(RhExp[0][0]       > MRh[0][0]) iAlarmMIN=1;
        if(RhExp[Ndata-1][0] < MRh[Nw-1][0]) iAlarmMAX=1;
    }
    else{
        istep=-1;
        Nflag0=Ndata-1;
        if(RhExp[Ndata-1][0] > MRh[0][0]) iAlarmMIN=1;
        if(RhExp[0][0]       < MRh[Nw-1][0]) iAlarmMAX=1;
    }
    if(iAlarmMIN!=0){
        msg="WLmin > "+QString::number(MRh[0][0]);
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
        return;
    }
    if(iAlarmMAX!=0){
        msg="WLmax < "+QString::number(MRh[Nw-1][0]);
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
        return;
    }
    // resampling RhExp like weight
    int jmax=1;
    if(ic==3 || ic==5) jmax=2;
    for(int j=1;j<=jmax;j++){
        Nflag=Nflag0;
        for(i=0;i<Nw;i++){
            WL=MRh[i][0];
            if(i>0)
                Wa=WL-(MRh[i][0]-MRh[i-1][0])/2.;
            else
                Wa=WL-(MRh[i+1][0]-MRh[i][0])/2.;
            if(i<Nw-1)
                Wb=WL+(MRh[i+1][0]-MRh[i][0])/2.;
            else
                Wb=WL+(MRh[i][0]-MRh[i-1][0])/2.;
            n=-1;
            do{
                if(RhExp[Nflag][0] >= Wa){
                    n++;
                    xy[0][n]=RhExp[Nflag][0];
                    xy[1][n]=RhExp[Nflag][j];
                }
                Nflag=Nflag+istep;
            } while(RhExp[Nflag][0] <= Wb && Nflag>=0 && Nflag<= Ndata-1);
            //      printf("Wa=%f Wb=%f n=%d \n",Wa,Wb,n);
            Nflag=Nflag-istep;
            Yfin=0.;
            if(n==-1){
                n=0;
                xy[0][n]=RhExp[Nflag][0];
                xy[1][n]=RhExp[Nflag][j];
                n=1;
                Nflag=Nflag+istep;
                xy[0][n]=RhExp[Nflag][0];
                xy[1][n]=RhExp[Nflag][j];
            }
            if(n==0){
                //	printf("Nflag=%d\n",Nflag);
                if(xy[0][n]<WL){
                    if(Nflag+istep >=0 && Nflag+istep <= Ndata-1){
                        //	    printf("caso 1\n");
                        Nflag=Nflag+istep;
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                    else{
                        //	    printf("caso 2\n");
                        n=0;
                        xy[0][n]=RhExp[Nflag-istep][0];
                        xy[1][n]=RhExp[Nflag-istep][j];
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                }
                else{
                    if(Nflag-istep >=0 && Nflag-istep <= Ndata-1){
                        //	    printf("caso 3\n");
                        n=0;
                        xy[0][n]=RhExp[Nflag-istep][0];
                        xy[1][n]=RhExp[Nflag-istep][j];
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                    else{
                        //	    printf("caso 4\n");
                        Nflag=Nflag+1;
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                }
            }
            if(n==1){
                Yfin=xy[1][0]+(xy[1][1]-xy[1][0])/(xy[0][1]-xy[0][0])*(WL-xy[0][0]);
            }
            if(n==2){
                DET=xy[0][0]*xy[0][0]*(xy[0][1]-xy[0][2])-xy[0][0]*(xy[0][1]*xy[0][1]-xy[0][2]*xy[0][2]);
                DET=DET+xy[0][1]*xy[0][1]*xy[0][2]-xy[0][2]*xy[0][2]*xy[0][1];
                A=xy[1][0]*(xy[0][1]-xy[0][2])-xy[0][0]*(xy[1][1]-xy[1][2]);
                A=(A+xy[1][1]*xy[0][2]-xy[1][2]*xy[0][1])/DET;
                B=xy[0][0]*xy[0][0]*(xy[1][1]-xy[1][2])-xy[1][0]*(xy[0][1]*xy[0][1]-xy[0][2]*xy[0][2]);
                B=(B+xy[0][1]*xy[0][1]*xy[1][2]-xy[0][2]*xy[0][2]*xy[1][1])/DET;
                C=xy[0][0]*xy[0][0]*(xy[0][1]*xy[1][2]-xy[0][2]*xy[1][1]);
                C=C-xy[0][0]*(xy[0][1]*xy[0][1]*xy[1][2]-xy[0][2]*xy[0][2]*xy[1][1]);
                C=(C+xy[1][0]*(xy[0][1]*xy[0][1]*xy[0][2]-xy[0][2]*xy[0][2]*xy[0][1]))/DET;
                Yfin=A*WL*WL+B*WL+C;
            }
            if(n>2){
                parabolicFit(n,A,B,C);
                Yfin=A*WL*WL+B*WL+C;
            }
            for(int jj=0;jj<=n;jj++)
                // 	printf("%f %f\n",xy[0][jj],xy[1][jj]);
                //      if(n>=2) printf("A=%f B=%f C=%f\n",A,B,C);
                if(ic==10)
                    MRh[i][ic]=Yfin*MRh[i][2];
                else
                    MRh[i][ic+j-1]=Yfin;
            //      printf("WL=%f MRh[%d][%d]=%f \n",MRh[i][0],i,ic,MRh[i][ic+j-1]);
            //      char ch;
            //       if(ic==3 || ic==5){
            // 	std::cout << "To continue press a key & enter: ";
            // 	std::cin  >> ch;
            //       }
        }
    }
}


void RsRpEMA(int iWL, double theta, double phiS, int iTISopt){
    if(iWL==10) printf("->RsRpEMA: iWL=%d theta=%f phiS=%f iTISopt=%d\n",iWL,theta,phiS,iTISopt);
    double R2s,R2p,alpha,ThetaOL,ratPhi,phiP,Ratio=1.;
    complex<double> theta1,theta2;
    complex<double> N0(1.0,0.0);//air
    complex<double> theta0(theta,0.0);
    complex<double> N1(MRh[iWL][3],-MRh[iWL][4]);//over-layer
    RsRp(N0,N1,theta0,R1s,R1p,theta1);
    ThetaOL=real(theta1);
    ratPhi=cos(theta)/cos(ThetaOL);
    phiP=ratPhi*phiS;
    complex<double> N2(MRh[iWL][5],-MRh[iWL][6]);//metallic-layer
    RsRp(N1,N2,theta1,R2s,R2p,theta2);
    if(iTISopt==0){
        double SigmaS=FunFit(phiS*1000.);
        double SigmaP=FunFit(phiP*1000.);
        double TISp=exp(-pow(4.*Pig*SigmaP*cos(ThetaOL)/MRh[iWL][0]*MRh[iWL][3],TISc));
        double TISs=exp(-pow(4.*Pig*SigmaS*cos(ThetaOL)/MRh[iWL][0]*MRh[iWL][3],2.));
        ratPhi=pow(ratPhi,2.);//pow(phiP/phiS,2.)
        Ratio=ratPhi*TISs+(1.-ratPhi)*(TISs+TISp)/2.;
        Ratio=ratPhi*TISs+(1.-ratPhi)*sqrt(TISs*TISp);
        if(iWL==10) printf("iWL=%d theta=%f phiS=%f phiP=%f SigmaS=%f SigmaP=%f TISs=%f TISp=%f ratPhi=%f Ratio=%f\n",
               iWL,theta,phiS,phiP,SigmaS,SigmaP,TISs,TISp,ratPhi,Ratio);
    }
    else if(iTISopt==1){
        double lam=MRh[iWL][0]/MRh[iWL][3];
        double argTIS=pow(4.*Pig*Pa[2]/lam,2.);//*cos(ThetaOL)
        double TIS=exp(-argTIS);
        double ERFP=erf(Pa[1]*Pig*phiP*cos(ThetaOL)/lam);
        double ERFS=erf(Pa[1]*Pig*phiS/lam);
        Ratio=TIS+(1.-TIS)*ERFP*ERFS;
//        double ERFStoSIN;
//        if(ThetaOL<0.01)
//            ERFStoSIN=2.*sqrt(Pig)*Pa[1]*phiS/lam;
//        else
//            ERFStoSIN=ERFS/sin(ThetaOL);
//        double argEQ19=2./pow(Pig,0.5)*exp(2)*pow(Pa[2]/Pa[1],2)/pow(lam,2)/cos(ThetaOL)*ERFStoSIN;
//        Ratio=Pa[1]*exp(-pow(Pa[1]*phiP*Pig*cos(ThetaOL)/lam,2))*phiP*lam*Pig*(3.*cos(ThetaOL)+5.*cos(3.*ThetaOL))+
//                sqrt(Pig)*(pow(lam,2)+3.*pow(Pa[1]*Pig,2)+(-5.*pow(lam,2)+4.*pow(Pa[1]*Pig,2))*cos(2.*ThetaOL)+
//                pow(Pa[1]*Pig,2)*cos(4.*ThetaOL))*ERFP;
//        Ratio=TIS+Pa[0]*argTIS*argEQ19*Ratio;
//        if(Ratio>1 || Ratio<TIS){
//            printf("Pa[0]=%f Pa[1]=%f Pa[2]=%f ThetaOL=%f lam=%f\n",Pa[0],Pa[1],Pa[2],ThetaOL,lam);
//            printf("phiS=%f phiP=%f ERFP=%f ERFS=%f argTIS=%f TIS=%f argEQ19=%E Ratio=%f\n",
//                   phiS,phiP,ERFP,ERFS,argTIS,TIS,argEQ19,Ratio);
//            if(Ratio<TIS)
//                Ratio=TIS;
//            if(Ratio>1.)
//                Ratio=1.;
//        }
        //Ratio=TIScomplete(phiS,phiP, MRh[iWL][0]/MRh[iWL][3], ThetaOL);
    }
    alpha=4.0*Pig/MRh[iWL][0]*imag(-sqrt(N1*N1-N0*N0*sin(theta0)*sin(theta0)));
    REMAs=R1s+pow(1.0-R1s,2.0)*R2s/(exp(2.*alpha*d_layer)-R1s*R2s);
    RnsEMAs=R1s+pow(1.0-R1s,2.0)*R2s*Ratio/(exp(2.*alpha*d_layer)-R1s*R2s*Ratio);
    REMAp=R1p+pow(1.0-R1p,2.0)*R2p/(exp(2.*alpha*d_layer)-R1p*R2p);
    RnsEMAp=R1p+pow(1.0-R1p,2.0)*R2p*Ratio/(exp(2.*alpha*d_layer)-R1p*R2p*Ratio);
    REMAu=0.5*(REMAs+REMAp);
    RnsEMAu=0.5*(RnsEMAs+RnsEMAp);
    if(iWL==10) printf("\tREMAs=%f REMAp=%f REMAu=%f\n",REMAs,REMAp,REMAu);
    if(iWL==10) printf("\tRnsEMAs=%f RnsEMAp=%f RnsEMAu=%f\n",RnsEMAs,RnsEMAp,RnsEMAu);
}



void RsRp(complex<double> N1, complex<double> N2, complex<double> theta1,
          double& Rs,double& Rp, complex<double>& theta2){
    theta2=asin(N1*sin(theta1)/N2);
    Rs=pow(abs((N1*cos(theta1)-N2*cos(theta2))/(N1*cos(theta1)+N2*cos(theta2))),2.);
    Rp=pow(abs((N2*cos(theta1)-N1*cos(theta2))/(N2*cos(theta1)+N1*cos(theta2))),2.);
    //  if(Rs<0. || Rp<0. || Rs>1.0 || Rp>1.0){
    //   cout << "N1= " << N1 << "\n";
    //   cout << "N2= " << N2 << "\n";
    //   cout << "cos(theta1)= " << cos(theta1) << "\n";
    //   cout << "cos(theta2)= " << cos(theta2) << "\n";
    //  }
}


double FunFit(double xp){
    double sigFun=Pa[0]+Pa[1]*exp(-xp*xp/Pa[2]/Pa[2])+Pa[3]*exp(-xp*xp/Pa[4]/Pa[4]);
    return(sigFun);
}


double TIScomplete(double phiS, double phiP,double lambdaOL, double ThetaOL){
    double TIS1=exp(-pow(4.*Pig*Pa[2]*cos(ThetaOL)/lambdaOL,2.));
    double TIS2=exp(-pow(4.*Pig*Pa[4]*cos(ThetaOL)/lambdaOL,2.));
    double TISco=TIS1*TIS2+(1.-TIS1*TIS2)*(
            Pa[0]*erf(Pa[1]/lambdaOL*phiP*cos(ThetaOL))*erf(Pa[1]/lambdaOL*phiS)
      +(1.-Pa[0])*erf(Pa[3]/lambdaOL*phiP*cos(ThetaOL))*erf(Pa[3]/lambdaOL*phiS));
    return(TISco);
}



void FindRoot(int ix,int i,int iWarning){
    double dt=0.1,dlim,tol,t,t0,d1,d0;
    //char ch;

    if(ix==4){//k over
        dt=abs(MRh[i][ix])/10.;
        if(dt < 1.0E-08) dt=1.0E-08;
        t=MRh[i][ix];
    }
    else if(ix==5 || ix==6){//n or k metal
        dt=0.01;
        t=MRh[i][ix];
    }
    else if(ix==56){//n&k metal
        dt=0.01;
        t=1.0;
        Nini=MRh[i][5];
        Kini=MRh[i][6];
    }
    else {
        return;
    }
    dlim=0.00001;
    tol=dt/1000.;

    d1=FT(ix,i,t);
    t0=t;
    d0=d1;
    int istop=0;
    int irif=0;
    int inul=0;
    int imi=0;
    int ilim=1000;
    if(iWarning==1)
      printf("\n FindRoot: i=%d t=%e dt=%e dlim=%e\n",i,t,dt,dlim);
    if(abs(d1) <= dlim) istop=1;
    while(istop == 0){
        t=t+dt;
        d1=FT(ix,i,t);
        if(iWarning==1){
           printf("t: %e -> %e dt=%e\td: %e -> %e DELTAd=%e\n",t0,t,dt,d0,d1,d1-d0);
//           cout << "To continue press a key & enter: ";
//           cin  >> ch;
        }
        if(abs(d1) <= dlim){
            istop=1;
            //    printf("  STOP per raggiunto limite precisione!\n");
        }
        else if(d0*d1 <.0 && abs(d1) >= (10.*abs(d0))){
            t=t0;
            d1=d0;
            dt=dt/10.;
            //    printf("dt decimato = %e\n",dt);
            if(abs(dt) <= dlim){
                istop=1;
                //      printf(" STOP per raggiunto limite precisione!\n");
            }
        }
        else if(d0*d1 <.0 && abs(d1) < (10.*abs(d0))) {// Brent
            double a=t;
            double b=t0;
            double c=b;
            double fa=d1;
            double fb=d0;
            double fc=d0;
            int itmax=1000;
            double eps=3.e-9;
            double xm,ss,pp,q,rr,e=0.,d=0.,tol1;
            int iter=0;
            while(istop == 0.){
                iter=iter+1;
                if((fb>0. && fc>0.) || (fb<0. && fc<0.)){
                    c=a; //Rename a, b, c and adjust bounding interval d.
                    fc=fa;
                    d=b-a;
                    e=d;
                }
                if(abs(fc) < abs(fb)){
                    a=b;
                    b=c;
                    c=a;
                    fa=fb;
                    fb=fc;
                    fc=fa;
                }
                //      printf("\nBrent: N_iterazione=%d \n",iter);
                //   "a,b,c=",a,b,c
                //   "fa,fb,fc=",fa,fb,fc
                tol1=2.*eps*abs(b)+0.5*tol; // Convergence check.
                xm=.5*(c-b);
                if(abs(xm)<=tol1 || fb==0. || iter>=itmax){
                    //	printf("eps=%e tol=%e tol1=%e \n xm=%e\n",eps,tol,tol1,xm);
                    t=b;
                    d1=fb;
                    istop=1;
                    //        if(iter<itmax)
                    //          printf("STOP per raggiunto limite di precisione\n");
                    //        else
                    //          printf("STOP per raggiunto N.max iterazioni\n");
                }
                if(istop==0){
                    if(abs(e)>=tol1 && abs(fa)>abs(fb)){
                        ss=fb/fa; //Attempt inverse quadratic interpolation.
                        if(a==c){
                            pp=2.*xm*ss;
                            q=1.-ss;
                        }
                        else{
                            q=fa/fc;
                            rr=fb/fc;
                            pp=ss*(2.*xm*q*(q-rr)-(b-a)*(rr-1.));
                            q=(q-1.)*(rr-1.)*(ss-1.);
                        }
                        if(pp>0.) q=-q; //Check whether in bounds.
                        pp=abs(pp);
                        if(2.*pp<min(3.*xm*q-abs(tol1*q),abs(e*q))){
                            e=d; //Accept interpolation.
                            d=pp/q;
                            //            printf("...interpolazione quadratica...\n");
                        }
                        else{
                            //            printf("...bisezione..\n");
                            d=xm; //Interpolation failed, use bisection.
                            e=d;
                        }
                    }
                    else{ //Bounds decreasing too slowly, use bisection.
                        //          printf("...bisezione per velocizzare..\n");
                        d=xm;
                        e=d;
                    }
                    a=b; //Move last best guess to a.
                    fa=fb;
                    if(abs(d)>tol1) //Evaluate new trial root.
                        b=b+d;
                    else{
                        if(xm>0)
                            b=b+tol1;
                        else
                            b=b-tol1;
                    }
                    t=b;
                    fb=FT(ix,i,t);
                    //        printf("'b=%e,fb=%e \n",b,fb);
                }
            }
        }
        else if(abs(d1)<abs(d0)) {
            //      printf("  -> miglioramento\n");
            t0=t;
            d0=d1;
            irif=0;
            imi=imi+1;
        }
        else if(abs(d1)>abs(d0) && d0*d1>0.){
            //      printf("  -> peggioramento\n");
            t=t0;
            dt=-dt;
            irif=irif+1;
            imi=0;
            if(irif>2){
                if(abs(d0-d1)>0.3*dlim){
                    dt=dt/2.;
                    //          printf("Incremento dimezzato!\n");
                }
                else
                    irif=ilim+10;
            }
            d1=d0;
        }
        else
            inul=inul+1;
        if(irif>ilim || imi>ilim || inul>ilim){
            //      printf("C0: irif o imi o inul > 1000! -> stop\n");
            istop=1;
        }
    }
    d1=FT(ix,i,t);
    if(iWarning==1)
      printf("Root -> t=%e d1=%e\n",t,d1);
    //   std::cout << "To continue press a key & enter: ";
    //   std::cin  >> ch;
}



double FT(int ix,int i,double t){
    double DeltaRh;
    if(ix>=0 && ix<56){
        MRh[i][ix]=t;
        RsRpEMA(i,thetaService,0.,-1);
        DeltaRh=REMAu-MRh[i][10];
    }
    else{
        MRh[i][5]=Nini*t;
        MRh[i][6]=Kini*t;
        RsRpEMA(i,thetaService,0.,-1);
        DeltaRh=REMAu-MRh[i][10];
    }
    return DeltaRh;
}


//*****************Soiled Mirror Model *************************
//**************************************************************

void SMEQ::getFileRhSoiled(){
    QString fileRhSoil = QFileDialog::getOpenFileName(
        this,
        "Choose a R-Hemispherical exp. file",
        pathroot+"/Spectra");
    if(fileRhSoil.isEmpty())
        return;
    lineEdit_RhSoiled -> setText(fileRhSoil.section('/',-1));
    spada(20,fileRhSoil);
    printf("Loaded %s\n",(fileRhSoil.toStdString()).c_str());
    double meanVal=MeanComputing(20);
    lineEdit_meanRhSoiled->setText(QString::number(meanVal,'f',4));

    //Plot Rh
    Gplot(1,0,0,0,"R-hemispherical spectrum",10);

    //Plot Rh_soiled
    for(int ii=0;ii< Nw;ii++){
        MRh[ii][24]=(MRh[ii][20]-MRh[ii][10])/MRh[ii][10];
    }
    Gplot(1,1,5,0,"R-hemispherical spectrum",20);
}

void SMEQ::getFileRnsSoiled(){
    QString fileRns = QFileDialog::getOpenFileName(
        this,
        "Choose a Rnear_specular exp. file",
        pathroot+"/Spectra");
    if(fileRns.isEmpty())
        return;
    lineEdit_RnsSoiled -> setText(fileRns.section('/',-1));
    spada(22,fileRns);
    printf("Loaded %s\n",(fileRns.toStdString()).c_str());
    double meanVal=MeanComputing(22);
    lineEdit_meanRnsSoiled->setText(QString::number(meanVal,'f',4));

    //Plot Rns_soiled
    for(int ii=0;ii< Nw;ii++){
        MRh[ii][26]=(MRh[ii][22]-MRh[ii][12])/MRh[ii][12];
    }
    Gplot(2,1,5,0,"R-nearspecular spectrum",22);
}

void SMEQ::PlotRns2(){
    PlotRns();
    Gplot(2,1,5,0,"R-nearspecular spectrum",22);
}

void SMEQ::plotNRloss(){

    //Plot NormRh_loss
    Gplot(3,0,0,0,"Norm-Reflectance loss",24);

    //Plot NormRns_loss
    Gplot(3,1,1,0,"Norm-Reflectance loss",26);

    //Plot computed NormRns_loss
    calcPlot();
}

void SMEQ::calcPlot(){
    if(iLock!=0)
        return;
    //model
    double S=doubleSpinBox_S->value();
    double sigmaSoil=doubleSpinBox_sigmaSoil->value();
    double C=0;//double C=doubleSpinBox_C->value();
    double D=doubleSpinBox_D->value();
    double k=doubleSpinBox_k->value();
    double L=doubleSpinBox_L->value();
    double B=doubleSpinBox_B->value();
    double E=doubleSpinBox_E->value();
    //double Ymin=doubleSpinBox_Ymin->value();
    //double Ymax=doubleSpinBox_Ymax->value();
    double WLmin=doubleSpinBox_Xmin->value();
    double WLmax=doubleSpinBox_Xmax->value();
    double WL,A;
    double stdDev=0.;
    int Ndafi=0;
    double thetaRns=dSB_thetaRns->value();
    thetaRns=thetaRns*Pig/180.;
    for(int i=0;i<Nw;i++){
        WL=MRh[i][0];
        A=MRh[i][24];
        MRh[i][27]=fmodel(WL,thetaRns,A,S,sigmaSoil,C,D,k,L,B,E);
        if(WL>=WLmin && WL<=WLmax){
            stdDev=stdDev+pow(MRh[i][27]-MRh[i][26],2.);
            Ndafi++;
        }
    }
    stdDev=sqrt(stdDev/Ndafi);
    lineEdit_stdDev->setText(QString::number(stdDev));
    Gplot(3,1,3,1,"Normalized reflectance loss",27);
}


void SMEQ::bestFit(){
    //best fit of normalized-reflectance loss
    iLock=1;
    double WLmin=doubleSpinBox_Xmin->value();
    double WLmax=doubleSpinBox_Xmax->value();
    int Ndafit=0;
    int Nmin=Nw;
    int Nmax=0;
    for(int i=0;i<Nw;i++){
        if(MRh[i][0]>=WLmin && MRh[i][0]<=WLmax){
            Ndafit++;
            Nmin=min(Nmin,i);
            Nmax=max(Nmax,i);
        }
    }
    printf("Ndafit= %d Nmin= %d Nmax= %d\n",Ndafit,Nmin,Nmax);
    Qt::CheckState state;
    for(int i=0; i<7; i++)
        ParFit[i][1]=0.;
    //fit parameters
    ParFit[0][0]=doubleSpinBox_S->value();
    ParFit[1][0]=doubleSpinBox_sigmaSoil->value();
    ParFit[2][0]=doubleSpinBox_D->value();
    ParFit[3][0]=doubleSpinBox_k->value();
    ParFit[4][0]=doubleSpinBox_L->value();
    ParFit[5][0]=doubleSpinBox_B->value();
    ParFit[6][0]=doubleSpinBox_E->value();
    thetaService=dSB_thetaRns->value();
    thetaService=thetaService*Pig/180.;//to be pass to fcn
    //ParFit[7][0]=doubleSpinBox_C->value();
    int n=0;
    state = checkBox_S -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[0][1]=1.;
    }
    state = checkBox_sigma -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[1][1]=1.;
    }
    state = checkBox_D -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[2][1]=1.;
    }
    state = checkBox_k -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[3][1]=1.;
    }
    state = checkBox_L -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[4][1]=1.;
    }
    state = checkBox_B -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[5][1]=1.;
    }
    state = checkBox_E -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[6][1]=1.;
    }
    int m=Ndafit;
    int lwa=m*n+5*n+m;
    vector <int> iwa(n);
    vector <double> x(n),fvec(m),wa(lwa);
    double tol=sqrt(dpmpar(1));
    int j=0;
    for(int i=0; i<7; i++){
        if(ParFit[i][1]>0.5){
            x[j]=ParFit[i][0];
            printf("Initial x[%d]= %f Nparameter= %d\n",j,x[j],i);
            j++;
        }
    }
    pTF2[0].NdatFit=Ndafit;
    pTF2[0].Nmin=Nmin;
    pTF2[0].Nmax=Nmax;
    int info=0;
    if(m>0 && n>0 && Ndafit>0){
        printf("Launch bestFit with %d parameters and %d data .... ",n,m);
        info=lmdif1(fcn, &pTF2, m, n, x.data(),fvec.data(), tol, iwa.data(), wa.data(), lwa);
        InfoFit(info);
    }
    else{
        printf("Nmis=0 || m=0 || n=0 then RETURN!");
        iLock=0;
        return;
    }
    printf(" done with info=%d\n",info);
    j=0;
    for(int i=0; i<7; i++){
        if(ParFit[i][1]>0.5){
            ParFit[i][0]=x[j];
            printf("Final x[%d]= %f Nparameter= %d\n",j,x[j],i);
            j++;
        }
    }
    doubleSpinBox_S    ->setValue(ParFit[0][0]);
    doubleSpinBox_sigmaSoil->setValue(ParFit[1][0]);
    doubleSpinBox_D    ->setValue(ParFit[2][0]);
    doubleSpinBox_k    ->setValue(ParFit[3][0]);
    doubleSpinBox_L    ->setValue(ParFit[4][0]);
    doubleSpinBox_B    ->setValue(ParFit[5][0]);
    doubleSpinBox_E    ->setValue(abs(ParFit[6][0]));
    iLock=0;
    calcPlot();
}



void SMEQ::RnsPrediction(){
    double thetaPre=dSB_thetaPrediction->value();
    thetaPre=thetaPre*Pig/180.;
    double phiS=dSB_Phi->value();
    phiS=phiS/1000.;
    TISc=doubleSpinBox_TIS->value();
    //soil model
    double S=doubleSpinBox_S->value();
    double sigmaSoil=doubleSpinBox_sigmaSoil->value();
    double C=0;//double C=doubleSpinBox_C->value();
    double D=doubleSpinBox_D->value();
    double k=doubleSpinBox_k->value();
    double L=doubleSpinBox_L->value();
    double B=doubleSpinBox_B->value();
    double E=doubleSpinBox_E->value();
    double sigmaClean=dSB_sigmaClean->value();
    double WL;
    Pa[0]=sigmaClean;//sigma (nm)
    printf("->RnsPrediction: thetaPre(rad)=%f phiS(rad)=%f sigmaSoil(nm)=%f\n",thetaPre,phiS,sigmaSoil);
    for(int i=0;i< Nw;i++){
        WL=MRh[i][0];
        MRh[i][27]=fmodel(WL,thetaPre,MRh[i][24],S,sigmaSoil,C,D,k,L,B,E);
        RsRpEMA(i, thetaPre, phiS, 0);//RnsCleanPrediction
        MRh[i][13]=RnsEMAu;
        MRh[i][28]=RnsEMAu*(1.+MRh[i][27]);//RnsSoiledPrediction
    }
    double meanVal=MeanComputing(28);
    lineEdit_meanRnsSoilPre->setText(QString::number(meanVal,'f',4));
    Gplot(3,1,6,1,"Normalized reflectance loss",27);
    Gplot(2,1,2,1,"R-nearspecular spectrum",28);
    Gplot(2,1,5,1,"R-nearspecular spectrum",13);
}


double fmodel(double wl,double theta,double A,double S,double sigma,double C,double D,double k,double L,double B,double E){
    double pi=3.14159265359;
    //double theRad=theta/180.*pi;
    double Le=L*sqrt(pow(cos(theta),2.)+pow(E*sin(theta),2.));
    double y=(A-S*(1-exp(-pow((4*pi*sigma*pow(cos(theta),C)/wl),2.)))-D/(1+exp(-2*(wl-Le)/k)))/pow(cos(theta),B);
    return(y);
}


int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag){
    //fit SoilMirModel
    struct pointToFit2 *pTF2 = (struct pointToFit2 *)p;
    int Nmin=pTF2[0].Nmin;
    int Nmax=pTF2[0].Nmax;
    double S,sigma,C,D,k,L,B,E;
    int j=0;
    if(ParFit[0][1]>0.5){
        S=x[j];
        j++;
    }
    else
        S=ParFit[0][0];
    if(ParFit[1][1]>0.5){
        sigma=x[j];
        j++;
    }
    else
        sigma=ParFit[1][0];
    if(ParFit[2][1]>0.5){
        D=x[j];
        j++;
    }
    else
        D=ParFit[2][0];
    if(ParFit[3][1]>0.5){
        k=x[j];
        j++;
    }
    else
        k=ParFit[3][0];
    if(ParFit[4][1]>0.5){
        L=x[j];
        j++;
    }
    else
        L=ParFit[4][0];
    if(ParFit[5][1]>0.5){
        B=x[j];
        j++;
    }
    else
        B=ParFit[5][0];
    if(ParFit[6][1]>0.5){
        E=x[j];
        j++;
    }
    else
        E=ParFit[6][0];
    C=0;
    int iv=0;
    double WL;
    for(j=Nmin;j<=Nmax;j++){
        WL=MRh[j][0];
        fvec[iv]=MRh[j][26]-fmodel(WL,thetaService,MRh[j][24],S,sigma,C,D,k,L,B,E);
        iv++;
    }
    return(0);
}
