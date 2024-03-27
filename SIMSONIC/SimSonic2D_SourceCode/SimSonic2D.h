/* Copyright (C) 2003-2012 BOSSY Emmanuel.

This file is part of the SimSonic suite.

SimSonic is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimSonic is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SimSonic.  If not, see <http://www.gnu.org/licenses/>. */

#define PREC_GRILLE float
#define PREC_CARTE unsigned char

/**************************************************************
 *************** Dimensions de la simulation ******************
 **************************************************************/

int X;
int Y;
int W;
int XPhotoAc;
int YPhotoAc;

/**************************************************************
 ************ Grilles de calcul du bloc central ***************
 **************************************************************/

PREC_GRILLE *Vx=NULL,*Vy=NULL,*Txx=NULL,*Tyy=NULL,*Txy=NULL;
PREC_GRILLE *Ux=NULL,*Uy=NULL;
PREC_GRILLE *MaxVx=NULL,*MaxVy=NULL,*MaxTxx=NULL,*MaxTyy=NULL,*MaxTxy=NULL;
PREC_GRILLE *MaxEc=NULL,*MaxV=NULL,*MaxCurl=NULL,*MaxDiv=NULL;
PREC_GRILLE *Vx_Old=NULL,*Vy_Old=NULL;

/**************************************************************
 ****************** Grilles des bords PML *********************
 **************************************************************/

PREC_GRILLE *Vx_Yp_X=NULL,*Vx_Yp_Y=NULL,*Vy_Yp_X=NULL,*Vy_Yp_Y=NULL;
PREC_GRILLE *Vx_Ym_X=NULL,*Vx_Ym_Y=NULL,*Vy_Ym_X=NULL,*Vy_Ym_Y=NULL;
PREC_GRILLE *Vx_Xp_X=NULL,*Vx_Xp_Y=NULL,*Vy_Xp_X=NULL,*Vy_Xp_Y=NULL;
PREC_GRILLE *Vx_Xm_X=NULL,*Vx_Xm_Y=NULL,*Vy_Xm_X=NULL,*Vy_Xm_Y=NULL;
PREC_GRILLE *Txx_Yp_X=NULL,*Txx_Yp_Y=NULL,*Tyy_Yp_X=NULL,*Tyy_Yp_Y=NULL,*Txy_Yp_X=NULL,*Txy_Yp_Y=NULL;
PREC_GRILLE *Txx_Ym_X=NULL,*Txx_Ym_Y=NULL,*Tyy_Ym_X=NULL,*Tyy_Ym_Y=NULL,*Txy_Ym_X=NULL,*Txy_Ym_Y=NULL;
PREC_GRILLE *Txx_Xp_X=NULL,*Txx_Xp_Y=NULL,*Tyy_Xp_X=NULL,*Tyy_Xp_Y=NULL,*Txy_Xp_X=NULL,*Txy_Xp_Y=NULL;
PREC_GRILLE *Txx_Xm_X=NULL,*Txx_Xm_Y=NULL,*Tyy_Xm_X=NULL,*Tyy_Xm_Y=NULL,*Txy_Xm_X=NULL,*Txy_Xm_Y=NULL;

/**************************************************************
 ****************** Grilles des coins PML *********************
 **************************************************************/

PREC_GRILLE *Vx_XmYp_X=NULL,*Vx_XmYp_Y=NULL,*Vy_XmYp_X=NULL,*Vy_XmYp_Y=NULL;
PREC_GRILLE *Vx_XmYm_X=NULL,*Vx_XmYm_Y=NULL,*Vy_XmYm_X=NULL,*Vy_XmYm_Y=NULL;
PREC_GRILLE *Vx_XpYp_X=NULL,*Vx_XpYp_Y=NULL,*Vy_XpYp_X=NULL,*Vy_XpYp_Y=NULL;
PREC_GRILLE *Vx_XpYm_X=NULL,*Vx_XpYm_Y=NULL,*Vy_XpYm_X=NULL,*Vy_XpYm_Y=NULL;
PREC_GRILLE *Txx_XmYp_X=NULL,*Txx_XmYp_Y=NULL,*Tyy_XmYp_X=NULL,*Tyy_XmYp_Y=NULL,*Txy_XmYp_X=NULL,*Txy_XmYp_Y=NULL;
PREC_GRILLE *Txx_XmYm_X=NULL,*Txx_XmYm_Y=NULL,*Tyy_XmYm_X=NULL,*Tyy_XmYm_Y=NULL,*Txy_XmYm_X=NULL,*Txy_XmYm_Y=NULL;
PREC_GRILLE *Txx_XpYp_X=NULL,*Txx_XpYp_Y=NULL,*Tyy_XpYp_X=NULL,*Tyy_XpYp_Y=NULL,*Txy_XpYp_X =NULL,*Txy_XpYp_Y=NULL;
PREC_GRILLE *Txx_XpYm_X=NULL,*Txx_XpYm_Y=NULL,*Tyy_XpYm_X=NULL,*Tyy_XpYm_Y=NULL,*Txy_XpYm_X=NULL,*Txy_XpYm_Y=NULL;

/**************************************************************
 ************** Grilles tampons des bords PML *****************
 **************************************************************/

PREC_GRILLE *Vi_BufferBordBlocI=NULL;
PREC_GRILLE *Tij_BufferBordBlocI=NULL;
PREC_GRILLE *Vi_BufferBordCoinJp=NULL;
PREC_GRILLE *Vi_BufferBordCoinJm=NULL;
PREC_GRILLE *Tjj_BufferBordCoinJp=NULL;
PREC_GRILLE *Tjj_BufferBordCoinJm=NULL;

/**************************************************************
 ************** Grilles tampons des coins PML *****************
 **************************************************************/

PREC_GRILLE *Vj_BufferCoinBordI=NULL;
PREC_GRILLE *Vi_BufferCoinBordJ=NULL;
PREC_GRILLE *Tij_BufferCoinBordI=NULL;
PREC_GRILLE *Tij_BufferCoinBordJ=NULL;

/**************************************************************
 *************** Tampons physiquement alloués *****************
 **************************************************************/

PREC_CARTE  *Carte_Tampon=NULL;
PREC_GRILLE *BufferBord01=NULL;
PREC_GRILLE *BufferBord02=NULL;
PREC_GRILLE *BufferBord03=NULL;
PREC_GRILLE *BufferBord04=NULL;

/**************************************************************
 ****************** Carte du bloc central *********************
 **************************************************************/

PREC_CARTE *Carte=NULL;
/* Ajouts de la version 2D_2006.03.16 */
PREC_CARTE *CartePhotoAc=NULL;
int ModePhotoAc=0;

/**************************************************************
 **************** Paramètres du bloc central ******************
 **************************************************************/

double *TableLegerete=NULL;
double *TableC11=NULL;
double *TableC22=NULL;
double *TableC12=NULL;
double *TableC33=NULL;
double *TableX11=NULL;
double *TableX22=NULL;
double *TableX12=NULL;
double *TableX33=NULL;
double *TableQModel=NULL;
double *AttenuationPML=NULL;
double *TableCii=NULL;
double *TableCjj=NULL;
double *TableCij=NULL;
double *TableCijij=NULL;
/* Ajout de la version 2D_2006.03.16 */
double *TableCoeffPhotoAc=NULL;

/**************************************************************
 ******************* Paramètres généraux **********************
 **************************************************************/

double H,Duree,Vmax,VmaxPML,DeltaT,CoeffDeltaT;
double DeltaTsurH;
int IndiceTempsMax,IndiceTemps;
double PeriodImage;
double PeriodRecepteurs;
int IndiceTempsImage;
int IndiceTempsRecepteurs;
int IndiceTempsAffichageTempsRestant;
int NumeroImage;
int TypeAbsorption;
double CoeffPML;
char Repertoire[500] ="";
char ImageFolder[500];
int FoundImageFolder=0;
int IsSkipFileCreated=0;
int IsDoneFileCreated=0;
int Is_TR_mode_with_initial_condition=0;
time_t TempsDebut,TempsCourant;
double TempsRestant;
int RcvDwnSmplng;
int SourceType;
int ConditionBordXm;
int ConditionBordXp;
int ConditionBordYp;
int ConditionBordYm;
int SnapshotV,SnapshotVx,SnapshotVy,SnapshotTxx,SnapshotTyy,SnapshotTxy,SnapshotEc;
int SnapshotMaxV,SnapshotMaxVx,SnapshotMaxVy,SnapshotMaxTxx,SnapshotMaxTyy,SnapshotMaxTxy,SnapshotMaxEc;
int SnapshotCurl,SnapshotDiv,SnapshotMaxCurl,SnapshotMaxDiv,SnapshotUx,SnapshotUy;
int NumberSelectedSnapshots=0;
int *SnapshotsTimeSelection=NULL;


/**************************************************************
 ***************** Gestion des traducteurs ********************
 **************************************************************/

double *SignalPhotoAc=NULL;

struct BARRETTE_EMISSION{
	char Normale;
	int XBar,YBar;
	int NbElements,Pitch,Largeur,Apodisation,Focale;
	double Direction;
	double VitesseMilieuEmission;
	int TypeSignal;
	char NomSignal[200];
	double Freq,T0,Bandwidth;
	double *Signal;
	/* Ajout de la version 2D_2005.10.04 */
	int DurationPt;
};

struct BARRETTE_RECEPTION{
	char Nom[200];
	char Normale;
	int XBar,YBar;
	int NbElements,Pitch,Largeur;
	int DurationPt;
	int ActiveMedium;
	double *Signaux;
};

int NbBarrettesEmissionTxx;
struct BARRETTE_EMISSION *EmissionTxx=NULL;
int NbBarrettesEmissionTyy;
struct BARRETTE_EMISSION *EmissionTyy=NULL;
int NbBarrettesEmissionTxy;
struct BARRETTE_EMISSION *EmissionTxy=NULL;
int NbBarrettesEmissionVx;
struct BARRETTE_EMISSION *EmissionVx=NULL;
int NbBarrettesEmissionVy;
struct BARRETTE_EMISSION *EmissionVy=NULL;
int NbBarrettesReceptionTxx;
struct BARRETTE_RECEPTION *ReceptionTxx=NULL;
int NbBarrettesReceptionTyy;
struct BARRETTE_RECEPTION *ReceptionTyy=NULL;
int NbBarrettesReceptionTxy;
struct BARRETTE_RECEPTION *ReceptionTxy=NULL;
int NbBarrettesReceptionVx;
struct BARRETTE_RECEPTION *ReceptionVx=NULL;
int NbBarrettesReceptionVy;
struct BARRETTE_RECEPTION *ReceptionVy=NULL;

int NbBarrettesSourceFilesTxx;
struct BARRETTE_RECEPTION *BarrettesSourceFilesTxx= NULL;
int NbBarrettesSourceFilesTyy;
struct BARRETTE_RECEPTION *BarrettesSourceFilesTyy= NULL;
int NbBarrettesSourceFilesTxy;
struct BARRETTE_RECEPTION *BarrettesSourceFilesTxy= NULL;
int NbBarrettesSourceFilesVx;
struct BARRETTE_RECEPTION *BarrettesSourceFilesVx= NULL;
int NbBarrettesSourceFilesVy;
struct BARRETTE_RECEPTION *BarrettesSourceFilesVy= NULL;

/**************************************************************
 **************** Prototypes des fonctions ********************
 **************************************************************/

/**************************************************************
 ******************* Calculs principaux ***********************
 **************************************************************/

void CalculVitesses(void);
void CalculContraintes(void);
void CalculViscosite(void);
void CalculBordsVitesses(void);
void CalculBordsContraintes(void);
void CalculViscositeBords(void);
void CalculViscositeBordsContraintes(void);
void CalculVitessesQModel(void);
void CalculContraintesQModel(void);
void CalculPMLBordContraintes(PREC_GRILLE *Vi_I,PREC_GRILLE *Vi_J,
	PREC_GRILLE *Vj_I,PREC_GRILLE *Vj_J,PREC_GRILLE *Tii_I,PREC_GRILLE *Tii_J,
	PREC_GRILLE *Tjj_I,PREC_GRILLE *Tjj_J,PREC_GRILLE *Tij_I,PREC_GRILLE *Tij_J,
	int J,const char *TypeBord);
void CalculPMLBordVitesses(PREC_GRILLE *Vi_I,PREC_GRILLE *Vi_J,
	PREC_GRILLE *Vj_I,PREC_GRILLE *Vj_J,PREC_GRILLE *Tii_I,PREC_GRILLE *Tii_J,
	PREC_GRILLE *Tjj_I,PREC_GRILLE *Tjj_J,PREC_GRILLE *Tij_I,PREC_GRILLE *Tij_J,
	int J,const char *TypeBord);
void CalculPMLCoinContraintes(PREC_GRILLE *Vi_I,PREC_GRILLE *Vi_J,
	PREC_GRILLE *Vj_I,PREC_GRILLE *Vj_J,PREC_GRILLE *Tii_I,PREC_GRILLE *Tii_J,
	PREC_GRILLE *Tjj_I,PREC_GRILLE *Tjj_J,PREC_GRILLE *Tij_I,PREC_GRILLE *Tij_J,
	const char *TypeCoin);
void CalculPMLCoinVitesses(PREC_GRILLE *Vi_I,PREC_GRILLE *Vi_J,
	PREC_GRILLE *Vj_I,PREC_GRILLE *Vj_J,PREC_GRILLE *Tii_I,PREC_GRILLE *Tii_J,
	PREC_GRILLE *Tjj_I,PREC_GRILLE *Tjj_J,PREC_GRILLE *Tij_I,PREC_GRILLE *Tij_J,
	const char *TypeCoin);

/**************************************************************
 ***************** Calculs intermédiaires *********************
 **************************************************************/

void SauvegardeVitesseOld(void);
void ComputeDisplacement(void);
void UpdateMaxSnapshot(void);

/**************************************************************
 ********************* Termes sources *************************
 **************************************************************/

void SourcesContraintes(void);
void SourcesVitesses(void);
void SourcesPhotoAc(void);

/**************************************************************
 ******************* Lecture d'entrées ************************
 **************************************************************/

void LectureParametres(void);
void LectureConditionsInitiales(char *NomFichierCondIni,PREC_GRILLE *GrilleLue,int DimX,int DimY);
void LectureSignal(struct BARRETTE_EMISSION *BarretteCourante,int NbPtsSignal);

/**************************************************************
 ****************** Ecriture de fichiers **********************
 **************************************************************/

void SauvegardeRecepteurs(void);
void SauvegardeSnapshot(int NumeroImageCourante);
void RecordMaxSnapshot(void);
void EnregistrementBarretteReception(PREC_GRILLE *GrilleCourante,struct BARRETTE_RECEPTION *BarretteCourante,int DimX,int DimY);

/**************************************************************
 ********************* Gestion mémoire ************************
 **************************************************************/

void AllocationGrilles(void);
void LiberationMemoire(void);

/**************************************************************
 *********************** Affichage ****************************
 **************************************************************/

void AffichageEndTime(void);
void AffichageStartTime(void);
void AffichageTempsRestant(void);
void AffichageDureeSimul(void);

/**************************************************************
 ***************** Manipulation de signaux ********************
 **************************************************************/

void GenerationSignal(int TypeSignal,double Freq,double T0,double Bandwidth,double *Signal);
double MaxAbsSignal(double* Signal,int NbptsSignal);

/**************************************************************
 ***************** Gestion des traducteurs ********************
 **************************************************************/

void LectureParametresBarretteEmission(FILE *Fich,char *LigneCourante,struct BARRETTE_EMISSION *BarretteCourante);
void LectureParametresBarretteReception(FILE *Fich,char *LigneCourante,struct BARRETTE_RECEPTION *BarretteCourante);
void LectureParametresBarretteSourceFile(FILE *Fich,char *LigneCourante,struct BARRETTE_RECEPTION *BarretteCourante);
int ValiditeBarretteEmission(struct BARRETTE_EMISSION *BarretteCourante,int LimX,int LimY);
int ValiditeBarretteReception(struct BARRETTE_RECEPTION *BarretteCourante,int LimX,int LimY);
void ActivationBarretteEmission(PREC_GRILLE *GrilleCourante,struct BARRETTE_EMISSION *BarretteCourante,int DimX,int DimY);
void ActivationBarretteSourceFile(PREC_GRILLE *GrilleCourante,struct BARRETTE_RECEPTION *BarretteCourante,int DimX,int DimY);
void MAJBarretteRecepteurs(void);


/**************************************************************
 ******************* Gestion des tampons **********************
 **************************************************************/

void MiseAJourContraintesTamponBordPML(const char *TypeBord);
void MiseAJourVitessesTamponBordPML(const char *TypeBord);
void MiseAJourContraintesTamponCoinPML(const char *TypeCoin);
void MiseAJourVitessesTamponCoinPML(const char *TypeCoin);
void ReMiseAZeroTampons(void);
void AdressageBuffer(void);
