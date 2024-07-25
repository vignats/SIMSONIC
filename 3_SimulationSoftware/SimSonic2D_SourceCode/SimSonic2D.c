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

/* Standard C libraries. */
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
/* OpenMP. */
#include <omp.h>
/* Functions prototypes and global variables declarations. */
#include "SimSonic2D.h"
/* Functions prototypes for PML computations. */
#include "PMLCoin2DVitesses.c"
#include "PMLCoin2DContraintes.c"
#include "PMLBordContraintes.c"
#include "PMLBordVitesses.c"

#define max(a, b)  (((a) > (b)) ? (a) : (b))

#define VERSION "2D.2012.04.26"

/**************************************************************
 ************************** Main ******************************
 **************************************************************/
int main(int argc,char *argv[])
{
    /* Various flags and initialization. */
    int FoundPhotoAcMapFileName=0;
    int FoundMapFileName=0;

    /* Generic loop index. */
    int LoopIdx=0;

    /* Strings for handling file names. */
    char MapFileName[500]="";
    char ParameterFileName[500]="";
    char InitialConditionsFileName[500]="";
    char PhotoacousticMapFileName[500]="";
    char VersionFileName[500]="";
    char SkipFileName[500]="";
    char DoneFileName[500]="";

    /* Strings used to read lines in text files. */
    char CurrentLine[500]="";

    /* Pointers for handling files. */
    FILE* ParameterFile=NULL;
    FILE* MapFile=NULL;
    FILE* PhotoacousticMapFile=NULL;
    FILE* InitialConditionsFile=NULL;
    FILE* VersionFile=NULL;
    FILE* SkipFile=NULL;
    FILE* DoneFile=NULL;

    /* Tests if program is called with at least one argument.
    If so, retrieves the working directory.
    If not, outputs a file with version info. */
    if(argc==1)
    {
        /* Write a flag file with version info in file name and file content. */
        strcpy(VersionFileName,VERSION);
        strcat(VersionFileName,""".ver");
        if((VersionFile=fopen(VersionFileName,"w+t"))!=NULL)
        {
            fprintf(VersionFile,"%s",VERSION);
            fclose(VersionFile);
        }
        LiberationMemoire();
        exit(1);
    }
    else
        strcpy(Repertoire,argv[1]);

    TempsDebut=time(NULL);


    /* Looks for a skip file. */
    strcpy(SkipFileName,Repertoire);
    strcat(SkipFileName,"""skip");
    if((SkipFile=fopen(SkipFileName,"rb"))!=NULL)
    {
        printf("skip file found...\n");
        exit(1);
    }

    /* Writes a flag file with version info in file name and file content. */
    strcpy(VersionFileName,Repertoire);
    strcat(VersionFileName,VERSION);
    strcat(VersionFileName,""".ver");
    if((VersionFile= fopen(VersionFileName,"w+t"))!=NULL)
    {
        fprintf(VersionFile,"%s",VERSION);
        fclose(VersionFile);
    }

    /* Allocates and initializes tables of physical properties */
    TableLegerete=(double*)calloc(256,sizeof(double));
    TableC11=(double*)calloc(256,sizeof(double));
    TableC22=(double*)calloc(256,sizeof(double));
    TableC12=(double*)calloc(256,sizeof(double));
    TableC33=(double*)calloc(256,sizeof(double));
    TableX11=(double*)calloc(256,sizeof(double));
    TableX22=(double*)calloc(256,sizeof(double));
    TableX12=(double*)calloc(256,sizeof(double));
    TableX33=(double*)calloc(256,sizeof(double));
    TableQModel=(double*)calloc(256,sizeof(double));
    TableCoeffPhotoAc=(double*)calloc(256,sizeof(double));

    for(LoopIdx=0; LoopIdx<256; LoopIdx++)
    {
        TableLegerete[LoopIdx]=1;
        TableC11[LoopIdx]=2.25;
        TableC22[LoopIdx]=2.25;
        TableC12[LoopIdx]=2.25;
        TableC33[LoopIdx]=0;
        TableX11[LoopIdx]=0;
        TableX22[LoopIdx]=0;
        TableX12[LoopIdx]=0;
        TableX33[LoopIdx]=0;
        TableQModel[LoopIdx]=0;
        TableCoeffPhotoAc[LoopIdx]=0;
    }

    /* Opens Parametres.ini2D to look for map file name. */
    strcpy(ParameterFileName,Repertoire);
    strcat(ParameterFileName,"""Parametres.ini2D");

    if((ParameterFile=fopen(ParameterFileName,"rb"))==NULL)
    {
        strcpy(ParameterFileName,Repertoire);
        strcat(ParameterFileName,"""Parameters.ini2D");

        if((ParameterFile=fopen(ParameterFileName,"rb"))==NULL)
        {

            printf("Parameters.ini2D NOT FOUND ...\n");
            exit(1);
        }
    }
    /* Looks for a map filename (optional for the user). */
    fseek(ParameterFile,0L,SEEK_SET);
    while(fgets(CurrentLine,500,ParameterFile)!=NULL)
        if(strstr(CurrentLine,"Carte map2D file name")!=NULL)
        {
            /* Loads the map's filename */
            sscanf(CurrentLine+31,"%s",MapFileName);

            if(!strcmp(MapFileName,""))
                FoundMapFileName=0;
            else
                FoundMapFileName=1;
        }

    /* Close Parametres.ini2D after initial search for map's filename. */
    fclose(ParameterFile);

    /* Use default filename "carte.map2D" if no map's filename was found. */
    if(!FoundMapFileName)
    {
        strcpy(MapFileName,Repertoire);
        strcat(MapFileName,"""carte.map2D");
    }

    /* Opens map file. */
    if((MapFile=fopen(MapFileName,"rb"))==NULL)
    {

        strcpy(MapFileName,Repertoire);
        strcat(MapFileName,"""Geometry.map2D");

        if((MapFile=fopen(MapFileName,"rb"))==NULL)
        {

            printf("Problem opening map file: File not found...\n");
            LiberationMemoire();
            exit(1);
        }
    }
    /* Reads map dimensions. */
    fread(&X,sizeof(int),1,MapFile);
    fread(&Y,sizeof(int),1,MapFile);

    /* Search in the parameter file for indication of a file containing a photoacoustic map.
    If found, checks for dimensions consistency with material map file. Finally, a flag ModePhotoAc
    is turned on or off. */
    ParameterFile=fopen(ParameterFileName,"rb");
    fseek(ParameterFile,0L,SEEK_SET);
    while(fgets(CurrentLine,500,ParameterFile)!=NULL)
        if(strstr(CurrentLine,"Photoacoustic Map")!=NULL)
        {
            sscanf(CurrentLine+31,"%s",PhotoacousticMapFileName);
            if(!strcmp(PhotoacousticMapFileName,""))
                FoundPhotoAcMapFileName=0;
            else
                FoundPhotoAcMapFileName=1;
        }
    fclose(ParameterFile);
    if(FoundPhotoAcMapFileName)
    {
        if((PhotoacousticMapFile=fopen(PhotoacousticMapFileName, "rb" ))==NULL)
        {
            printf("Problem with the photoacoustic map:file not found...\n");
            LiberationMemoire();
            exit(1);
        }

        fread(&XPhotoAc,sizeof(int),1,PhotoacousticMapFile);
        fread(&YPhotoAc,sizeof(int),1,PhotoacousticMapFile);

        fclose(ParameterFile);
        if (!((XPhotoAc==X)&&(YPhotoAc==Y)))
        {
            printf("Error: sizes of the photoacoustic map and medium map do not match!\n");
            LiberationMemoire();
            exit(1);
        }
        ModePhotoAc=1;
    }

    /* Reads ParameterFile. (nb: uses ModePhotoAc flag). */
    LectureParametres();

    /* Dynamically allocates grids, all dimensions now being known. */
    AllocationGrilles();
    AdressageBuffer();

    /* Reads data from mapfile, and closes mapfile. */
    fread(Carte,sizeof(PREC_CARTE),X*Y,MapFile);
    fclose(MapFile);

    /* If necesssary, reads data from photoacoustic mapfile. */
    if (ModePhotoAc)
    {
        PhotoacousticMapFile=fopen(PhotoacousticMapFileName, "rb" );
        fread(&XPhotoAc,sizeof(int),1,PhotoacousticMapFile);
        fread(&YPhotoAc,sizeof(int),1,PhotoacousticMapFile);
        fread(CartePhotoAc,sizeof(PREC_CARTE),X*Y,PhotoacousticMapFile);
        fclose( PhotoacousticMapFile );
    }

    /* Looks for initial conditions, and reads if they exist. */
    strcpy(InitialConditionsFileName,Repertoire);
    strcat(InitialConditionsFileName,"""Vx0.snp2D");
    if(((InitialConditionsFile=fopen(InitialConditionsFileName,"rb"))!=NULL))
    {
        LectureConditionsInitiales(InitialConditionsFileName,Vx,X+1,Y);
        fclose(InitialConditionsFile);
    }
    strcpy(InitialConditionsFileName,Repertoire);
    strcat(InitialConditionsFileName,"""Vy0.snp2D");
    if(((InitialConditionsFile =fopen(InitialConditionsFileName,"rb"))!=NULL))
    {
        LectureConditionsInitiales(InitialConditionsFileName,Vy,X,Y+1);
        fclose(InitialConditionsFile);
    }
    strcpy(InitialConditionsFileName,Repertoire);
    strcat(InitialConditionsFileName,"""Txx0.snp2D");
    if(((InitialConditionsFile =fopen(InitialConditionsFileName,"rb"))!=NULL))
    {
        LectureConditionsInitiales(InitialConditionsFileName,Txx,X,Y);
        fclose(InitialConditionsFile);
    }
    strcpy(InitialConditionsFileName,Repertoire);
    strcat(InitialConditionsFileName,"""T220.snp2D");
    if(((InitialConditionsFile =fopen(InitialConditionsFileName,"rb"))!=NULL))
    {
        LectureConditionsInitiales(InitialConditionsFileName,Tyy,X,Y);
        fclose(InitialConditionsFile);
    }
    strcpy(InitialConditionsFileName,Repertoire);
    strcat(InitialConditionsFileName,"""T120.snp2D");
    if(((InitialConditionsFile =fopen(InitialConditionsFileName,"rb"))!=NULL))
    {
        LectureConditionsInitiales(InitialConditionsFileName,Txy,X+1,Y+1);
        fclose(InitialConditionsFile);
    }

    /* Defines "attenuation" profile in PMLs. */
#pragma omp parallel for schedule(guided)
    for(LoopIdx=0; LoopIdx<(2*W+1); LoopIdx++)
    {
        AttenuationPML[LoopIdx]=8.63E-2*((double)CoeffPML)*VmaxPML/(((double)W)*H)*((double)(LoopIdx*LoopIdx))/((double)(4*W*W));
        /* 8.63E-2 corresponds to 1/2*3/2*ln(10)/20 ... c.f. p.416 INRIA report). */
    }

    /* Writes a skip file. */
    if(IsSkipFileCreated)
        if((SkipFile=fopen(SkipFileName,"w+t"))!=NULL)
        {
            fclose(SkipFile);
        }

    /**************************************************************
     ************************ Main loop ***************************
     **************************************************************/


    NumeroImage=0;

    if(Is_TR_mode_with_initial_condition)
    {
        printf("Time Reversal Mode \n");
#pragma omp parallel for schedule(guided)
        for(LoopIdx=0; LoopIdx<(X+1)*Y; LoopIdx++)
            Vx[LoopIdx]=-Vx[LoopIdx];
#pragma omp parallel for schedule(guided)
        for(LoopIdx=0; LoopIdx<X*(Y+1); LoopIdx++)
            Vy[LoopIdx]=-Vy[LoopIdx];
        /* Bloc central. */
        CalculContraintes();
        /* Bords. */
        CalculBordsContraintes();
    }

    printf("Running %s\n",Repertoire);
    AffichageStartTime();

    for(IndiceTemps=0; IndiceTemps<IndiceTempsMax; IndiceTemps++)
    {
        /*********************************************/
        /*Calculation of displacement velocity fields*/
        /*********************************************/

        /* If a viscous model is used for absorption, this saves current velocity fields before updating. */
        if(TypeAbsorption==1)
            SauvegardeVitesseOld();

        CalculVitesses(); /* main grid */
        CalculBordsVitesses(); /* boundaries */

        /*add source terms. Needs to be placed AFTER calculations,to work with forced values,
        and to avoid Q-model to attenuate the source...*/
        SourcesVitesses();


        /*********************************************/
        /********* Calculation of stress fields ******/
        /*********************************************/


        CalculContraintes();/* main grid */
        CalculBordsContraintes();/* boundaries */


        /*if viscous model is used for absorption*/
        if(TypeAbsorption==1)
        {
            CalculViscosite();
            /*CalculViscositeBords(); TO DO....  */
        }


        /*If photoacoustic sources are used*/
        if(ModePhotoAc==1)
            SourcesPhotoAc();

        SourcesContraintes();

        /*********************************************/
        /********* Calculation of optional fields ******/
        /*********************************************/

        if(SnapshotUx||SnapshotUy)
            ComputeDisplacement();

        /*For Max Snapshot*/
        if(SnapshotMaxTxx||SnapshotMaxTyy||SnapshotMaxTxy||SnapshotMaxVx||SnapshotMaxVy||
                SnapshotMaxV||SnapshotMaxEc||SnapshotMaxDiv||SnapshotMaxCurl)
            UpdateMaxSnapshot();

        /*********************************************/
        /****** Files recording and various display **/
        /*********************************************/

        if(!((IndiceTemps)%RcvDwnSmplng))
        {
            MAJBarretteRecepteurs();
        }

        /*Records snapshots*/


        if(NumberSelectedSnapshots)
        {
            if(IndiceTemps==SnapshotsTimeSelection[NumeroImage])
            {
                NumeroImage++;
                SauvegardeSnapshot(NumeroImage);
                RecordMaxSnapshot();
            }
        }
        else
            if(!((IndiceTemps + 1)%IndiceTempsImage))
            {
                NumeroImage++;
                SauvegardeSnapshot(NumeroImage);
                RecordMaxSnapshot();
            }

        // ajouter une condition permettant les images selectionnes...


        /*Records signal on receivers*/
        if(!((IndiceTemps + 1)%IndiceTempsRecepteurs))
        {
            SauvegardeRecepteurs();
        }

        /* if necessary, displays time info to console*/
        if(argc > 2)
            if(!((IndiceTemps)%IndiceTempsAffichageTempsRestant))
                AffichageTempsRestant();

        if(argc==4)
        {
            if(IsSkipFileCreated)
                remove(SkipFileName);
            LiberationMemoire();
            exit(1);
        }
    }

    RecordMaxSnapshot();
    AffichageEndTime();
    AffichageDureeSimul();

    SauvegardeRecepteurs();
    LiberationMemoire();

    /* write a "done" file */
    strcpy(DoneFileName,Repertoire);
    strcat(DoneFileName,"""done");
    if(IsDoneFileCreated)
        if((DoneFile=fopen(DoneFileName,"w+t")) !=NULL)
        {
            fclose(DoneFile);
        }

    return 0;
}


/************************************************************************************************/
/**********************************   Fonctions  ************************************************/
/************************************************************************************************/

/************************************** CALCUL **************************************************/

/* Calculation of DISPLACEMENT VELOCITY in the MAIN GRID */
void CalculVitesses(void)
{
    /* Indice (x,y). */
    int i,j,index=0;
    /* Paramètre effectif de légereté. */
    double Legerete;
    double FacteurQModel;

#pragma omp parallel for private(j,Legerete,index) schedule(guided)
    for(i=1; i<X; i++)
    {
        for(j=0; j<Y; j++)
        {
            index=i*Y+j;
            Legerete=2./(1./TableLegerete[Carte[index-Y]]+1./TableLegerete[Carte[index]]);
            Vx[index]+=DeltaTsurH*Legerete*((Txx[index]-Txx[index-Y])+(Txy[index+i+1]-Txy[index+i]));
        }
    }
#pragma omp parallel for private(j,Legerete,index) schedule(guided)
    for(i=0; i<X; i++)
    {
        for(j=1; j<Y; j++)
        {
            index=i*Y+j;
            Legerete=2./(1./TableLegerete[Carte[index-1]]+1./TableLegerete[Carte[index]]);
            Vy[index+i]+=DeltaTsurH*Legerete*((Txy[index+i+Y+1]-Txy[index+i])+(Tyy[index]-Tyy[index-1]));
        }
    }

    /**************************************************************************/
    /*********************** Attenuating Factor (QModel) **********************/
    /**************************************************************************/

    if(TypeAbsorption==2)
    {
#pragma omp parallel for private(j,FacteurQModel) schedule(guided)
        for(i=1; i<X; i++)
        {
            for(j=0; j<Y; j++)
            {
                FacteurQModel=sqrt(TableQModel[Carte[(i-1)*Y+j]]*TableQModel[Carte[i*Y+j]]);
                Vx[i*Y+j] *=FacteurQModel;
            }
        }
#pragma omp parallel for private(j,FacteurQModel) schedule(guided)
        for(i=0; i<X; i++)
        {
            for(j=1; j<Y; j++)
            {
                FacteurQModel=sqrt(TableQModel[Carte[i*Y+j-1]]*TableQModel[Carte[i*Y+j]]);
                Vy[i*(Y+1)+j] *=FacteurQModel;
            }
        }
    }
}

/* Calculation of STRESS in the main grid */
void CalculContraintes(void)
{
    int i,j,index=0,position=0;
    double C33;
    double FacteurQModel;
    PREC_GRILLE DeltaVx,DeltaVy;

#pragma omp parallel for private(j,index,DeltaVx,DeltaVy,position) schedule(guided)
    for(i=0; i<X; i++)
    {
        for(j=0; j<Y; j++)
        {
            index=i*Y+j;
            DeltaVx=Vx[index+Y]-Vx[index];
            DeltaVy=Vy[index+i+1]-Vy[index+i];
            position=Carte[index];
            Txx[index]+=DeltaTsurH*(TableC11[position]*DeltaVx+TableC12[position]*DeltaVy);
            Tyy[index]+=DeltaTsurH*(TableC12[position]*DeltaVx+TableC22[position]*DeltaVy);
        }
    }
#pragma omp parallel for private(j,C33,index) schedule(guided)
    for(i=1; i<X; i++)
    {
        for(j=1; j<Y; j++)
        {
            index=i*Y+j;
            if(TableC33[Carte[index-Y-1]]*
                    TableC33[Carte[index-Y]]*
                    TableC33[Carte[index-1]]*
                    TableC33[Carte[index]])
            {
                C33=4./(1./TableC33[Carte[index-Y-1]]+
                        1./TableC33[Carte[index-Y]]+
                        1./TableC33[Carte[index-1]]+
                        1./TableC33[Carte[index]]);
            }
            else
            {
                C33=0.0;
            }
            Txy[index+i]+=DeltaTsurH*C33*((Vx[index]-Vx[index-1])+(Vy[index+i]-Vy[index+i-Y-1]));
        }
    }

    /**************************************************************************/
    /*********************** Attenuating Factor QModel ************************/
    /**************************************************************************/

    if(TypeAbsorption==2)
    {
#pragma omp parallel for private(j,index) schedule(guided)
        for(i=0; i<X; i++)
        {
            for(j=0; j<Y; j++)
            {
                index=i*Y+j;
                Txx[index]*=TableQModel[Carte[index]];
                Tyy[index]*=TableQModel[Carte[index]];
            }
        }
#pragma omp parallel for private(j,FacteurQModel,index) schedule(guided)
        for(i=1; i<X; i++)
        {
            for(j=1; j<Y; j++)
            {
                index=i*Y+j;
                FacteurQModel=sqrt(sqrt(TableQModel[Carte[index-Y-1]]*
                                        TableQModel[Carte[index-Y]]*
                                        TableQModel[Carte[index-1]]*
                                        TableQModel[Carte[index]]));
                Txy[index+i]*=FacteurQModel;
            }
        }
    }
}

/* Calculation of VISCOUS STRESS TERM in the MAIN GRID */
void CalculViscosite(void)
{
    int i,j;
    double X33;
    double DerivH=1./H;

#pragma omp parallel for private(j) schedule(guided)
    for(i=0; i<X; i++)
    {
        for(j=0; j < Y; j++)
        {
            Txx[i*Y+j] +=DerivH*(
                             TableX11[Carte[i*Y+j]]*(
                                 (Vx[(i+1)*Y+j]-Vx[(i)*Y+j])-(Vx_Old[(i+1)*Y+j]-Vx_Old[(i)*Y+j]))+
                             TableX12[Carte[i*Y+j]]*(
                                 (Vy[i*(Y+1)+j+1]-Vy[i*(Y+1)+j])-(Vy_Old[i*(Y+1)+j+1]-Vy_Old[i*(Y+1)+j])));
            Tyy[i*Y+j] +=DerivH*(
                             TableX12[Carte[i*Y+j]]*(
                                 (Vx[(i+1)*Y+j]-Vx[(i)*Y+j])-(Vx_Old[(i+1)*Y+j]-Vx_Old[(i)*Y+j]))+
                             TableX22[Carte[i*Y+j]]*(
                                 (Vy[i*(Y+1)+j+1]-Vy[i*(Y+1)+j])-(Vy_Old[i*(Y+1)+j+1]-Vy_Old[i*(Y+1)+j])));
        }
    }

#pragma omp parallel for private(j,X33) schedule(guided)
    for(i=1; i<X; i++)
    {
        for(j=1; j<Y; j++)
        {
            if(TableX33[Carte[(i-1)*Y+j-1]]*
                    TableX33[Carte[(i-1)*Y+j]]*
                    TableX33[Carte[i*Y+j-1]]*
                    TableX33[Carte[i*Y+j]])
            {
                X33=4./(1./TableX33[Carte[(i-1)*Y+j-1]]+
                        1./TableX33[Carte[(i-1)*Y+j]]+
                        1./TableX33[Carte[i*Y+j-1]]+
                        1./TableX33[Carte[i*Y+j]]);
            }
            else
            {
                X33=0.0;
            }
            Txy[i*(Y+1)+j] +=DerivH*X33*(
                                 (Vx[i*Y+j]-Vx[i*Y+j-1])+(Vy[(i)*(Y+1)+j]-Vy[(i-1)*(Y+1)+j]) -
                                 ((Vx_Old[i*Y+j]-Vx_Old[i*Y+j-1])+(Vy_Old[(i)*(Y+1)+j]-Vy_Old[(i-1)*(Y+1)+j]))
                             );
        }
    }

}
void CalculBordsVitesses(void)
{
    int j; /*Indice (x,y)*/
    double Legerete;  /*Paramètre effectif de légèreté*/
    double FacteurQModel;
    /* (0 : PML, 1 : Symetrie, 2 : Libre, 3 : Rigide, 4: anti-symétrie) */
    /* Bord Xp. */
    switch(ConditionBordXp)
    {
    case 0:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<Y; j++)
        {
            Legerete=TableLegerete[Carte[(X-1)*Y+j]];
            Vx[X*Y+j]+=DeltaTsurH*Legerete*(((Txx_Xp_X[j]+Txx_Xp_Y[j])-Txx[(X-1)*Y+j])+(Txy[X*(Y+1)+j+1]-Txy[X*(Y+1)+j]));
        }
        break;
    case 2:
    case 4:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<Y; j++)
        {
            Legerete=TableLegerete[Carte[(X-1)*Y+j]];
            Vx[X*Y+j]+=DeltaTsurH*Legerete*((0-2*Txx[(X-1)*Y+j])+(Txy[X*(Y+1)+j+1]-Txy[X*(Y+1)+j]));
        }
        break;
        /* Nul implicitement...*/
    case 1:
    case 3:
    default:
        break;
    }
    /* Bord Xm. */
    switch(ConditionBordXm)
    {
    case 0:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<Y; j++)
        {
            Legerete=TableLegerete[Carte[(0)*Y+j]];
            Vx[0*Y+j]+=DeltaTsurH*Legerete*(-((Txx_Xm_X[j]+Txx_Xm_Y[j])-Txx[0*Y+j])+(Txy[0*(Y+1)+j+1]-Txy[0*(Y+1)+j]));
        }
        break;
    case 2:
    case 4:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<Y; j++)
        {
            Legerete=TableLegerete[Carte[(0)*Y+j]];
            Vx[0*Y+j]+=DeltaTsurH*Legerete*(-(0-2*Txx[0*Y+j])+(Txy[0*(Y+1)+j+1]-Txy[0*(Y+1)+j]));
        }
        break;
        /* Nul implicitement ... */
    case 1:
    case 3:
    default:
        break;
    }
    /* Bord Yp. */
    switch(ConditionBordYp)
    {
    case 0:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<X; j++)
        {
            Legerete=TableLegerete[Carte[(j)*(Y)+Y-1]];
            Vy[j*(Y+1)+Y]+=DeltaTsurH*Legerete*((Txy[(j+1)*(Y+1)+Y]-Txy[(j)*(Y+1)+Y])+((Tyy_Yp_X[j]+Tyy_Yp_Y[j])-Tyy[j*(Y)+Y-1]));
        }
        break;
    case 2:
    case 4:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<X; j++)
        {
            Legerete=TableLegerete[Carte[(j)*(Y)+Y-1]];
            Vy[j*(Y+1)+Y]+=DeltaTsurH*Legerete*((Txy[(j+1)*(Y+1)+Y]-Txy[(j)*(Y+1)+Y])+(0-2*Tyy[j*(Y)+Y-1]));
        }
        break;
        /* Nul implicitement ... */
    case 1:
    case 3:
    default:
        break;
    }
    /* Bord Ym. */
    switch(ConditionBordYm)
    {
    case 0:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<X; j++)
        {
            Legerete=TableLegerete[Carte[(j)*(Y)+0]];
            Vy[j*(Y+1)+0]+=DeltaTsurH*Legerete*((Txy[(j+1)*(Y+1)+0]-Txy[(j)*(Y+1)+0])-((Tyy_Ym_X[j]+Tyy_Ym_Y[j])-Tyy[j*(Y)+0]));
        }
        break;
    case 2:
    case 4:
#pragma omp parallel for private(Legerete) schedule(guided)
        for(j=0; j<X; j++)
        {
            Legerete=TableLegerete[Carte[(j)*(Y)+0]];
            Vy[j*(Y+1)+0]+=DeltaTsurH*Legerete*((Txy[(j+1)*(Y+1)+0]-Txy[(j)*(Y+1)+0])-(0-2*Tyy[j*(Y)+0]));
        }
        break;
        /* Nul implicitement ...*/
    case 1:
    case 3:
    default:
        break;
    }

    /**************************************************************************/
    /*********************** Attenuating Factor QModel ************************/
    /**************************************************************************/
    if(TypeAbsorption==2)
    {
#pragma omp parallel for private(FacteurQModel) schedule(guided)
        for(j=0; j < Y; j++)
        {
            FacteurQModel=TableQModel[Carte[(X-1)*Y+j]];
            Vx[X*Y+j] *=FacteurQModel;
            FacteurQModel=TableQModel[Carte[(0)*Y+j]];
            Vx[0*Y+j] *=FacteurQModel;
        }
#pragma omp parallel for private(FacteurQModel) schedule(guided)
        for(j=0; j < X; j++)
        {
            FacteurQModel=TableQModel[Carte[(j)*Y+Y-1]];
            Vy[j*(Y+1)+Y] *=FacteurQModel;
            FacteurQModel=TableQModel[Carte[(j)*Y+0]];
            Vy[j*(Y+1)+0] *=FacteurQModel;
        }

    }
    /**************************************************************************/
    /*********************** Calcul des Bords PML *****************************/
    /**************************************************************************/
    if(ConditionBordXp==0)
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponBordPML("Xp");
        CalculPMLBordVitesses(Vx_Xp_X,Vx_Xp_Y,
                              Vy_Xp_X,Vy_Xp_Y,Txx_Xp_X,Txx_Xp_Y,
                              Tyy_Xp_X,Tyy_Xp_Y,Txy_Xp_X,Txy_Xp_Y,Y,"Xp");
    }
    if(ConditionBordYp==0)
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponBordPML("Yp");
        CalculPMLBordVitesses(Vy_Yp_Y,Vy_Yp_X,
                              Vx_Yp_Y,Vx_Yp_X,Tyy_Yp_Y,Tyy_Yp_X,
                              Txx_Yp_Y,Txx_Yp_X,Txy_Yp_Y,Txy_Yp_X,X,"Yp");
    }
    if(ConditionBordXm==0)
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponBordPML("Xm");
        CalculPMLBordVitesses(Vx_Xm_X,Vx_Xm_Y,
                              Vy_Xm_X,Vy_Xm_Y,Txx_Xm_X,Txx_Xm_Y,
                              Tyy_Xm_X,Tyy_Xm_Y,Txy_Xm_X,Txy_Xm_Y,Y,"Xm");
    }
    if(ConditionBordYm==0)
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponBordPML("Ym");
        CalculPMLBordVitesses(Vy_Ym_Y,Vy_Ym_X,
                              Vx_Ym_Y,Vx_Ym_X,Tyy_Ym_Y,Tyy_Ym_X,
                              Txx_Ym_Y,Txx_Ym_X,Txy_Ym_Y,Txy_Ym_X,X,"Ym");
    }
    /**************************************************************************/
    /*********************** Calcul des Coins PML *****************************/
    /**************************************************************************/
    if((ConditionBordXp==0) && (ConditionBordYp==0))
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponCoinPML("XpYp");
        CalculPMLCoinVitesses(Vx_XpYp_X,Vx_XpYp_Y,
                              Vy_XpYp_X,Vy_XpYp_Y,Txx_XpYp_X,Txx_XpYp_Y,
                              Tyy_XpYp_X,Tyy_XpYp_Y,Txy_XpYp_X,Txy_XpYp_Y,"XpYp");
    }
    if((ConditionBordXp==0) && (ConditionBordYm==0))
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponCoinPML("XpYm");
        CalculPMLCoinVitesses(Vx_XpYm_X,Vx_XpYm_Y,
                              Vy_XpYm_X,Vy_XpYm_Y,Txx_XpYm_X,Txx_XpYm_Y,
                              Tyy_XpYm_X,Tyy_XpYm_Y,Txy_XpYm_X,Txy_XpYm_Y,"XpYm");
    }
    if((ConditionBordXm==0) && (ConditionBordYp==0))
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponCoinPML("XmYp");
        CalculPMLCoinVitesses(Vx_XmYp_X,Vx_XmYp_Y,
                              Vy_XmYp_X,Vy_XmYp_Y,Txx_XmYp_X,Txx_XmYp_Y,
                              Tyy_XmYp_X,Tyy_XmYp_Y,Txy_XmYp_X,Txy_XmYp_Y,"XmYp");
    }
    if((ConditionBordXm==0) && (ConditionBordYm==0))
    {
        ReMiseAZeroTampons();
        MiseAJourContraintesTamponCoinPML("XmYm");
        CalculPMLCoinVitesses(Vx_XmYm_X,Vx_XmYm_Y,
                              Vy_XmYm_X,Vy_XmYm_Y,Txx_XmYm_X,Txx_XmYm_Y,
                              Tyy_XmYm_X,Tyy_XmYm_Y,Txy_XmYm_X,Txy_XmYm_Y,"XmYm");
    }
}

/* ToDo */
void CalculViscositeBordsContraintes(void)
{
}

void CalculBordsContraintes(void)
{
    int j;
    double C33;
    double FacteurQModel;
    /* (0 : PML, 1 : Symétrie, 2 : Libre, 3 : Rigide, 4: anti-symétrie) */
    /* Bord Xp. */
    switch(ConditionBordXp)
    {
    case 0:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<Y; j++)
        {
            if(TableC33[Carte[(X-1)*Y+j-1]]*TableC33[Carte[(X-1)*Y+j]])
                C33=2./(1./TableC33[Carte[(X-1)*Y+j-1]]+1./TableC33[Carte[(X-1)*Y+j]]);
            else
                C33=0.0;
            Txy[X*(Y+1)+j]+=DeltaTsurH*C33*((Vx[X*Y+j]-Vx[X*Y+j-1])+((Vy_Xp_X[j]+Vy_Xp_Y[j])-Vy[(X-1)*(Y+1)+j]));
        }
        break;
        /* Nul implicitement ... */
    case 1:
    case 2:
    default:
        break;
    case 3:
    case 4:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<Y; j++)
        {
            if(TableC33[Carte[(X-1)*Y+j-1]]*TableC33[Carte[(X-1)*Y+j]])
                C33=2./(1./TableC33[Carte[(X-1)*Y+j-1]]+1./TableC33[Carte[(X-1)*Y+j]]);
            else
                C33=0.0;
            Txy[X*(Y+1)+j]+=DeltaTsurH*C33*((Vx[X*Y+j]-Vx[X*Y+j-1])+(0-2*Vy[(X-1)*(Y+1)+j]));
        }
        break;
    }
    /* Bord Xm. */
    switch(ConditionBordXm)
    {
    case 0:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<Y; j++)
        {
            if(TableC33[Carte[(0)*Y+j-1]]*TableC33[Carte[(0)*Y+j]])
                C33=2./(1./TableC33[Carte[(0)*Y+j-1]]+1./TableC33[Carte[(0)*Y+j]]);
            else
                C33=0.0;
            Txy[0*(Y+1)+j]+=DeltaTsurH*C33*((Vx[0*Y+j]-Vx[0*Y+j-1])-((Vy_Xm_X[j]+Vy_Xm_Y[j])-Vy[(0)*(Y+1)+j]));
        }
        break;
        /* Nul implicitement ... */
    case 1:
    case 2:
    default:
        break;
    case 3:
    case 4:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<Y; j++)
        {
            if(TableC33[Carte[(0)*Y+j-1]]*TableC33[Carte[(0)*Y+j]])
                C33=2./(1./TableC33[Carte[(0)*Y+j-1]]+1./TableC33[Carte[(0)*Y+j]]);
            else
                C33=0.0;
            Txy[0*(Y+1)+j]+=DeltaTsurH*C33*((Vx[0*Y+j]-Vx[0*Y+j-1])-(0-2*Vy[(0)*(Y+1)+j]));
        }
        break;
    }
    /* Bord Yp. */
    switch(ConditionBordYp)
    {
    case 0:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<X; j++)
        {
            if(TableC33[Carte[(j-1)*Y+Y-1]]*TableC33[Carte[(j)*Y+Y-1]])
                C33=2./(1./TableC33[Carte[(j-1)*Y+Y-1]]+1./TableC33[Carte[(j)*Y+Y-1]]);
            else
                C33=0.0;
            Txy[j*(Y+1)+Y]+=DeltaTsurH*C33*(((Vx_Yp_X[j]+Vx_Yp_Y[j])-Vx[j*Y+Y-1])+(Vy[(j)*(Y+1)+Y]-Vy[(j-1)*(Y+1)+Y]));
        }
        break;
        /* Nul implicitement ... */
    case 1:
    case 2:
    default:
        break;
    case 3:
    case 4:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<X; j++)
        {
            if(TableC33[Carte[(j-1)*Y+Y-1]]*TableC33[Carte[(j)*Y+Y-1]])
                C33=2./(1./TableC33[Carte[(j-1)*Y+Y-1]]+1./TableC33[Carte[(j)*Y+Y-1]]);
            else
                C33=0.0;
            Txy[j*(Y+1)+Y]+=DeltaTsurH*C33*((0-2*Vx[j*Y+Y-1])+(Vy[(j)*(Y+1)+Y]-Vy[(j-1)*(Y+1)+Y]));
        }
        break;
    }
    /* Bord Ym/ */
    switch(ConditionBordYm)
    {
    case 0:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<X; j++)
        {
            if(TableC33[Carte[(j-1)*Y+0]]*TableC33[Carte[(j)*Y+0]])
                C33=2./(1./TableC33[Carte[(j-1)*Y+0]]+1./TableC33[Carte[(j)*Y+0]]);
            else
                C33=0.0;
            Txy[j*(Y+1)+0]+=DeltaTsurH*C33*(-((Vx_Ym_X[j]+Vx_Ym_Y[j])-Vx[j*Y+0])+(Vy[(j)*(Y+1)+0]-Vy[(j-1)*(Y+1)+0]));
        }
        break;
        /* Nul implicitement ... */
    case 1:
    case 2:
    default:
        break;
    case 3:
    case 4:
#pragma omp parallel for private(C33) schedule(guided)
        for(j=1; j<X; j++)
        {
            if(TableC33[Carte[(j-1)*Y+0]]*TableC33[Carte[(j)*Y+0]])
                C33=2./(1./TableC33[Carte[(j-1)*Y+0]]+1./TableC33[Carte[(j)*Y+0]]);
            else
                C33=0.0;
            Txy[j*(Y+1)+0] +=DeltaTsurH*C33*(-(0 - 2*Vx[j*Y+0]) + (Vy[(j)*(Y+1)+0] - Vy[(j-1)*(Y+1)+0]));
        }
        break;
    }
    /******************** Coin XpYp *******************/
    if((ConditionBordXp==0) && (ConditionBordYp==0))
    {
        C33=TableC33[Carte[(X-1)*Y+Y-1]];
        Txy[X*(Y+1)+Y] +=DeltaTsurH*C33*(((Vx_Yp_X[X]+Vx_Yp_Y[X]) - Vx[X*Y+Y-1])+
                                         ((Vy_Xp_X[Y]+Vy_Xp_Y[Y]) - Vy[(X-1)*(Y+1)+Y]));
    }
    if((ConditionBordXp==0) && ((ConditionBordYp==3)||(ConditionBordYp==4)))
    {
        C33=TableC33[Carte[(X-1)*Y+Y-1]];
        Txy[X*(Y+1)+Y] +=DeltaTsurH*C33*((0 - 2*Vx[X*Y+Y-1]) + ((Vy_Xp_X[Y]+Vy_Xp_Y[Y]) - Vy[(X-1)*(Y+1)+Y]));
    }
    if(((ConditionBordXp==3) || (ConditionBordXp==4)) && (ConditionBordYp==0))
    {
        C33=TableC33[Carte[(X-1)*Y+Y-1]];
        Txy[X*(Y+1)+Y] +=DeltaTsurH*C33*(((Vx_Yp_X[X]+Vx_Yp_Y[X]) - Vx[X*Y+Y-1]) + (0 - 2*Vy[(X-1)*(Y+1)+Y]));
    }
    if(((ConditionBordXp==3)||(ConditionBordXp==4)) && ((ConditionBordYp==3)||(ConditionBordYp==4)))
    {
        C33=TableC33[Carte[(X-1)*Y+Y-1]];
        Txy[X*(Y+1)+Y] +=DeltaTsurH*C33*((0 - 2*Vx[X*Y+Y-1]) + (0 - 2*Vy[(X-1)*(Y+1)+Y]));
    }
    /******************** Coin XpYm *******************/
    if((ConditionBordXp==0) && (ConditionBordYm==0))
    {
        C33=TableC33[Carte[(X-1)*Y+0]];
        Txy[X*(Y+1)+0] +=DeltaTsurH*C33*(-((Vx_Ym_X[X]+Vx_Ym_Y[X]) - Vx[X*Y+0]) + ((Vy_Xp_X[0]+Vy_Xp_Y[0]) - Vy[(X-1)*(Y+1)+0]));
    }
    if((ConditionBordXp==0) && ((ConditionBordYm==3)||(ConditionBordYm==4)))
    {
        C33=TableC33[Carte[(X-1)*Y+0]];
        Txy[X*(Y+1)+0] +=DeltaTsurH*C33*(-(0 - 2*Vx[X*Y+0]) + ((Vy_Xp_X[0]+Vy_Xp_Y[0]) - Vy[(X-1)*(Y+1)+0]));
    }
    if(((ConditionBordXp==3)||(ConditionBordXp==4)) && (ConditionBordYm==0))
    {
        C33=TableC33[Carte[(X-1)*Y+0]];
        Txy[X*(Y+1)+0] +=DeltaTsurH*C33*(-((Vx_Ym_X[X]+Vx_Ym_Y[X]) - Vx[X*Y+0]) + (0 - 2*Vy[(X-1)*(Y+1)+0]));
    }
    if(((ConditionBordXp==3)||(ConditionBordXp==4)) && ((ConditionBordYm==3)||(ConditionBordYm==4)))
    {
        C33=TableC33[Carte[(X-1)*Y+0]];
        Txy[X*(Y+1)+0] +=DeltaTsurH*C33*(-(0 - 2*Vx[X*Y+0]) + (0 - 2*Vy[(X-1)*(Y+1)+0]));
    }
    /******************** Coin XmYp *******************/
    if((ConditionBordXm==0) && (ConditionBordYp==0))
    {
        C33=TableC33[Carte[(0)*Y+Y-1]];
        Txy[0*(Y+1)+Y] +=DeltaTsurH*C33*(((Vx_Yp_X[0]+Vx_Yp_Y[0]) - Vx[0*Y+Y-1]) - ((Vy_Xm_X[Y]+Vy_Xm_Y[Y]) - Vy[(0)*(Y+1)+Y]));
    }
    if((ConditionBordXm==0) && ((ConditionBordYp==3)||(ConditionBordYp==4)))
    {
        C33=TableC33[Carte[(0)*Y+Y-1]];
        Txy[0*(Y+1)+Y] +=DeltaTsurH*C33*((0 - 2*Vx[0*Y+Y-1]) - ((Vy_Xm_X[Y]+Vy_Xm_Y[Y]) - Vy[(0)*(Y+1)+Y]));
    }
    if(((ConditionBordXm==3)||(ConditionBordXm==4)) && (ConditionBordYp==0))
    {
        C33=TableC33[Carte[(0)*Y+Y-1]];
        Txy[0*(Y+1)+Y] +=DeltaTsurH*C33*(((Vx_Yp_X[0]+Vx_Yp_Y[0]) - Vx[0*Y+Y-1]) - (0 - 2*Vy[(0)*(Y+1)+Y]));
    }
    if(((ConditionBordXm==3)||(ConditionBordXm==4)) && ((ConditionBordYp==3)||(ConditionBordYp==4)))
    {
        C33=TableC33[Carte[(0)*Y+Y-1]];
        Txy[0*(Y+1)+Y] +=DeltaTsurH*C33*((0 - 2*Vx[0*Y+Y-1]) - (0 - 2*Vy[(0)*(Y+1)+Y]));
    }
    /******************** Coin XmYm *******************/
    if((ConditionBordXm==0) && (ConditionBordYm==0))
    {
        C33=TableC33[Carte[(0)*Y+0]];
        Txy[0*(Y+1)+0] +=DeltaTsurH*C33*(-((Vx_Ym_X[0]+Vx_Ym_Y[0]) - Vx[0*Y+0]) - ((Vy_Xm_X[0]+Vy_Xm_Y[0]) - Vy[(0)*(Y+1)+0]));
    }
    if((ConditionBordXm==0) && ((ConditionBordYm==3)||(ConditionBordYm==4)))
    {
        C33=TableC33[Carte[(0)*Y+0]];
        Txy[0*(Y+1)+0] +=DeltaTsurH*C33*(-(0 - 2*Vx[0*Y+0]) - ((Vy_Xm_X[0]+Vy_Xm_Y[0]) - Vy[(0)*(Y+1)+0]));
    }
    if(((ConditionBordXm==3)||(ConditionBordXm==4)) && (ConditionBordYm==0))
    {
        C33=TableC33[Carte[(0)*Y+0]];
        Txy[0*(Y+1)+0] +=DeltaTsurH*C33*(-((Vx_Ym_X[0]+Vx_Ym_Y[0]) - Vx[0*Y+0]) - (0 - 2*Vy[(0)*(Y+1)+0]));
    }
    if(((ConditionBordXm==3)||(ConditionBordXm==4)) && ((ConditionBordYm==3)||(ConditionBordYm==4)))
    {
        C33=TableC33[Carte[(0)*Y+0]];
        Txy[0*(Y+1)+0] +=DeltaTsurH*C33*(-(0 - 2*Vx[0*Y+0]) - (0 - 2*Vy[(0)*(Y+1)+0]));
    }
    /**************************************************************************/
    /*********************** Attenuating Factor QModel ************************/
    /**************************************************************************/
    if(TypeAbsorption==2)
    {
#pragma omp parallel for private(FacteurQModel) schedule(guided)
        for(j=1; j < Y; j++)
        {
            FacteurQModel=sqrt(TableQModel[Carte[(X-1)*Y+j-1]]*TableQModel[Carte[(X-1)*Y+j]]);
            Txy[X*(Y+1)+j] *=FacteurQModel;
            FacteurQModel=sqrt(TableQModel[Carte[(0)*Y+j-1]]*TableQModel[Carte[(0)*Y+j]]);
            Txy[0*(Y+1)+j] *=FacteurQModel;
        }
#pragma omp parallel for private(FacteurQModel) schedule(guided)
        for(j=1; j < X; j++)
        {
            FacteurQModel=sqrt(TableQModel[Carte[(j-1)*Y+Y-1]]*TableQModel[Carte[(j)*Y+Y-1]]);
            Txy[j*(Y+1)+Y] *=FacteurQModel;
            FacteurQModel=sqrt(TableQModel[Carte[(j-1)*Y+0]]*TableQModel[Carte[(j)*Y+0]]);
            Txy[j*(Y+1)+0] *=FacteurQModel;
        }

        Txy[0*(Y+1)+0] *=TableQModel[Carte[(0)*Y+0]];
        Txy[0*(Y+1)+Y] *=TableQModel[Carte[(0)*Y+Y-1]];
        Txy[X*(Y+1)+0] *=TableQModel[Carte[(X-1)*Y+0]];
        Txy[X*(Y+1)+Y] *=TableQModel[Carte[(X-1)*Y+Y-1]];
    }
    /**************************************************************************/
    /*********************** Calcul des Bords PML *****************************/
    /**************************************************************************/
    if(ConditionBordXp==0)
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponBordPML("Xp");
        CalculPMLBordContraintes(Vx_Xp_X,Vx_Xp_Y,
                                 Vy_Xp_X,Vy_Xp_Y,Txx_Xp_X,Txx_Xp_Y,
                                 Tyy_Xp_X,Tyy_Xp_Y,Txy_Xp_X,Txy_Xp_Y,Y,"Xp");
    }

    if(ConditionBordYp==0)
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponBordPML("Yp");
        CalculPMLBordContraintes(Vy_Yp_Y,Vy_Yp_X,
                                 Vx_Yp_Y,Vx_Yp_X,Tyy_Yp_Y,Tyy_Yp_X,
                                 Txx_Yp_Y,Txx_Yp_X,Txy_Yp_Y,Txy_Yp_X,X,"Yp");
    }
    if(ConditionBordXm==0)
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponBordPML("Xm");
        CalculPMLBordContraintes(Vx_Xm_X,Vx_Xm_Y,
                                 Vy_Xm_X,Vy_Xm_Y,Txx_Xm_X,Txx_Xm_Y,
                                 Tyy_Xm_X,Tyy_Xm_Y,Txy_Xm_X,Txy_Xm_Y,Y,"Xm");
    }
    if(ConditionBordYm==0)
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponBordPML("Ym");
        CalculPMLBordContraintes(Vy_Ym_Y,Vy_Ym_X,
                                 Vx_Ym_Y,Vx_Ym_X,Tyy_Ym_Y,Tyy_Ym_X,
                                 Txx_Ym_Y,Txx_Ym_X,Txy_Ym_Y,Txy_Ym_X,X,"Ym");
    }
    /**************************************************************************/
    /*********************** Calcul des Coins PML *****************************/
    /**************************************************************************/
    if((ConditionBordXp==0) && (ConditionBordYp==0))
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponCoinPML("XpYp");
        CalculPMLCoinContraintes(Vx_XpYp_X,Vx_XpYp_Y,
                                 Vy_XpYp_X,Vy_XpYp_Y,Txx_XpYp_X,Txx_XpYp_Y,
                                 Tyy_XpYp_X,Tyy_XpYp_Y,Txy_XpYp_X,Txy_XpYp_Y,"XpYp");
    }
    if((ConditionBordXp==0) && (ConditionBordYm==0))
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponCoinPML("XpYm");
        CalculPMLCoinContraintes(Vx_XpYm_X,Vx_XpYm_Y,
                                 Vy_XpYm_X,Vy_XpYm_Y,Txx_XpYm_X,Txx_XpYm_Y,
                                 Tyy_XpYm_X,Tyy_XpYm_Y,Txy_XpYm_X,Txy_XpYm_Y,"XpYm");
    }
    if((ConditionBordXm==0) && (ConditionBordYp==0))
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponCoinPML("XmYp");
        CalculPMLCoinContraintes(Vx_XmYp_X,Vx_XmYp_Y,
                                 Vy_XmYp_X,Vy_XmYp_Y,Txx_XmYp_X,Txx_XmYp_Y,
                                 Tyy_XmYp_X,Tyy_XmYp_Y,Txy_XmYp_X,Txy_XmYp_Y,"XmYp");
    }
    if((ConditionBordXm==0) && (ConditionBordYm==0))
    {
        ReMiseAZeroTampons();
        MiseAJourVitessesTamponCoinPML("XmYm");
        CalculPMLCoinContraintes(Vx_XmYm_X,Vx_XmYm_Y,
                                 Vy_XmYm_X,Vy_XmYm_Y,Txx_XmYm_X,Txx_XmYm_Y,
                                 Tyy_XmYm_X,Tyy_XmYm_Y,Txy_XmYm_X,Txy_XmYm_Y,"XmYm");
    }
}
/**************************************** CALCULS INTERMEDIAIRES *************************************/

void SauvegardeVitesseOld(void)
{
    int i,j;
#pragma omp parallel for private(j) schedule(guided)
    for(i=0; i < (X+1); i++)
    {
        for(j=0; j < Y; j++)
        {
            Vx_Old[i*Y+j]=Vx[i*Y+j];
        }
    }
#pragma omp parallel for private(j) schedule(guided)
    for(i=0; i < X; i++)
    {
        for(j=0; j < (Y+1); j++)
        {
            Vy_Old[i*(Y+1)+j]=Vy[i*(Y+1)+j];
        }
    }

    /* A FAIRE: Sauvez les vitesses dans les PML activées....*/
}
void ComputeDisplacement(void) /* computes displacement step by step by integrating velocity */
{


    int DimX,DimY;
    int i;


    if(SnapshotUx)
    {

        DimX=X +1;
        DimY=Y;
#pragma omp parallel for schedule(guided)
        for(i=0; i < DimX*DimY; i++)
            Ux[i]+=DeltaT*Vx[i];
    }
    if(SnapshotUy)
    {

        DimX=X;
        DimY=Y + 1;
#pragma omp parallel for schedule(guided)
        for(i=0; i < DimX*DimY; i++)
            Uy[i]+=DeltaT*Vy[i];
    }

}
void UpdateMaxSnapshot(void)
{
    int DimX,DimY;
    int i,j;
    float Temp;
    double Legerete;
    double DerivH=1./H;

    if(SnapshotMaxTxx)
    {
        DimX=X;
        DimY=Y;
#pragma omp parallel for schedule(guided)
        for(i=0; i<DimX*DimY; i++)
        {
            if(fabs(Txx[i])>MaxTxx[i])
                MaxTxx[i]=fabs(Txx[i]);
        }
    }

    if(SnapshotMaxTyy)
    {
        DimX=X;
        DimY=Y;
#pragma omp parallel for schedule(guided)
        for(i=0; i< DimX*DimY; i++)
        {
            if(fabs(Tyy[i])>MaxTyy[i])
                MaxTyy[i]=fabs(Tyy[i]);
        }
    }

    if(SnapshotMaxTxy)
    {
        DimX=X+1;
        DimY=Y+1;
#pragma omp parallel for schedule(guided)
        for(i=0; i<DimX*DimY; i++)
        {
            if(fabs(Txy[i])>MaxTxy[i])
                MaxTxy[i]=fabs(Txy[i]);
        }
    }

    if(SnapshotMaxVx)
    {

        DimX=X+1;
        DimY=Y;
#pragma omp parallel for schedule(guided)
        for(i=0; i<DimX*DimY; i++)
        {
            if(fabs(Vx[i])>MaxVx[i])
                MaxVx[i]=fabs(Vx[i]);
        }
    }

    if(SnapshotMaxVy)
    {
        DimX=X;
        DimY=Y+1;
#pragma omp parallel for schedule(guided)
        for(i=0; i<DimX*DimY; i++)
        {
            if(fabs(Vy[i])>MaxVy[i])
                MaxVy[i]=fabs(Vy[i]);
        }
    }

    if(SnapshotMaxDiv)
    {
        DimX=X;
        DimY=Y;
#pragma omp parallel for private(j,Temp) schedule(guided)
        for(i=0; i<DimX; i++)
        {
            for(j=0; j<DimY; j++)
            {
                Temp=DerivH*((Vx[(i+1)*Y+j]-Vx[(i)*Y+j])+(Vy[i*(Y+1)+j+1]-Vy[i*(Y+1)+j]));
                if(fabs(Temp)>MaxDiv[i*DimY+j])
                    MaxDiv[i*DimY+j]=fabs(Temp);
            }
        }
    }

    if(SnapshotMaxCurl)
    {

        DimX=X-1;
        DimY=Y-1;
#pragma omp parallel for private(j,Temp) schedule(guided)
        for(i=0; i < DimX; i++)
        {
            for(j=0; j < DimY; j++)
            {

                Temp=DerivH*((Vy[(i+1)*(Y+1)+j+1]-Vy[i*(Y+1)+j+1])-(Vx[(i+1)*Y+j+1]-Vx[(i+1)*Y+j]));

                if (fabs(Temp)>MaxCurl[i*DimY+j])
                    MaxCurl[i*DimY+j]=fabs(Temp);
            }
        }
    }

    if(SnapshotMaxV)
    {

        DimX=X;
        DimY=Y;
#pragma omp parallel for private(j,Temp) schedule(guided)
        for(i=0; i < DimX; i++)
        {
            for(j=0; j < DimY; j++)
            {
                Temp=(float)sqrt(Vy[i*(Y+1)+j]*Vy[i*(Y+1)+j]+Vx[i*Y+j]*Vx[i*Y+j]);
                if (fabs(Temp)>MaxV[i*DimY+j])
                    MaxV[i*DimY+j]=fabs(Temp);
            }
        }

    }
    if(SnapshotMaxEc)
    {

        DimX=X;
        DimY=Y;
#pragma omp parallel for private(j,Temp,Legerete) schedule(guided)
        for(i=0; i < DimX; i++)
        {
            for(j=0; j < DimY; j++)
            {
                Legerete=(double)(TableLegerete[Carte[i*Y+j]]);
                Temp=(float)(0.5*(Vy[i*(Y+1)+j]*Vy[i*(Y+1)+j]+Vx[i*Y+j]*Vx[i*Y+j])/Legerete);
                if (fabs(Temp)>MaxEc[i*DimY+j])
                    MaxEc[i*DimY+j]=fabs(Temp);
            }
        }
    }

}


/************************ TERMES SOURCES **************************************/

void SourcesContraintes(void)
{
    int NumeroBarrette;
    struct BARRETTE_EMISSION *BarretteCourante;
    struct BARRETTE_RECEPTION *BarretteCouranteSourceFile;

    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesSourceFilesTxx; NumeroBarrette++)
    {
        BarretteCouranteSourceFile=&BarrettesSourceFilesTxx[NumeroBarrette];
        ActivationBarretteSourceFile(Txx,BarretteCouranteSourceFile,X,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesSourceFilesTyy; NumeroBarrette++)
    {
        BarretteCouranteSourceFile=&BarrettesSourceFilesTyy[NumeroBarrette];
        ActivationBarretteSourceFile(Tyy,BarretteCouranteSourceFile,X,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesSourceFilesTxy; NumeroBarrette++)
    {
        BarretteCouranteSourceFile=&BarrettesSourceFilesTxy[NumeroBarrette];
        ActivationBarretteSourceFile(Txy,BarretteCouranteSourceFile,X+1,Y+1);
    }


    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesEmissionTxx; NumeroBarrette++)
    {
        BarretteCourante=&EmissionTxx[NumeroBarrette];
        ActivationBarretteEmission(Txx,BarretteCourante,X,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesEmissionTyy; NumeroBarrette++)
    {
        BarretteCourante=&EmissionTyy[NumeroBarrette];
        ActivationBarretteEmission(Tyy,BarretteCourante,X,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesEmissionTxy; NumeroBarrette++)
    {
        BarretteCourante=&EmissionTxy[NumeroBarrette];
        ActivationBarretteEmission(Txy,BarretteCourante,X+1,Y+1);
    }
}
void SourcesVitesses(void)
{
    int NumeroBarrette;
    struct BARRETTE_EMISSION *BarretteCourante;
    struct BARRETTE_RECEPTION *BarretteCouranteSourceFile;

    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesSourceFilesVx; NumeroBarrette++)
    {
        BarretteCouranteSourceFile=&BarrettesSourceFilesVx[NumeroBarrette];
        ActivationBarretteSourceFile(Vx,BarretteCouranteSourceFile,X+1,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesSourceFilesVy; NumeroBarrette++)
    {
        BarretteCouranteSourceFile=&BarrettesSourceFilesVy[NumeroBarrette];
        ActivationBarretteSourceFile(Vy,BarretteCouranteSourceFile,X,Y+1);
    }


    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesEmissionVx; NumeroBarrette++)
    {
        BarretteCourante=&EmissionVx[NumeroBarrette];
        ActivationBarretteEmission(Vx,BarretteCourante,X+1,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesEmissionVy; NumeroBarrette++)
    {
        BarretteCourante=&EmissionVy[NumeroBarrette];
        ActivationBarretteEmission(Vy,BarretteCourante,X,Y+1);
    }
}

void SourcesPhotoAc(void)
{
    int i,j;
#pragma omp parallel for private(j) schedule(guided)
    for (i=0; i < X; i++)
    {
        for (j=0; j < Y; j++)
        {
            Txx[i*Y+j] -=TableCoeffPhotoAc[CartePhotoAc[i*Y+j]]*SignalPhotoAc[IndiceTemps]*DeltaT;
            Tyy[i*Y+j] -=TableCoeffPhotoAc[CartePhotoAc[i*Y+j]]*SignalPhotoAc[IndiceTemps]*DeltaT;
        }
    }

}

/**************   Gestion des Transducteurs ********************************/

void LectureParametresBarretteEmission(FILE *Fich, char *LigneCourante,struct BARRETTE_EMISSION *BarretteCourante)
{
    BarretteCourante->Signal=(double *) calloc (IndiceTempsMax,sizeof(double));
    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%i",&BarretteCourante->TypeSignal);
    if(BarretteCourante->TypeSignal)
    {
        if((BarretteCourante->TypeSignal>0))
        {
            sscanf(LigneCourante,"%i %lf %lf %lf",&BarretteCourante->TypeSignal,
                   &BarretteCourante->Freq,&BarretteCourante->T0,&BarretteCourante->Bandwidth);
            GenerationSignal(BarretteCourante->TypeSignal,BarretteCourante->Freq,
                             BarretteCourante->T0,BarretteCourante->Bandwidth,BarretteCourante->Signal);
        }
        else
        {
            /* add been added to keep compatiblity while extending folder management. Would not do it the
            	same way from scratch....*/
            char TempFileName[500]="";

            sscanf(LigneCourante,"%i %s",&BarretteCourante->TypeSignal,BarretteCourante->NomSignal);

            strcpy(TempFileName,BarretteCourante->NomSignal);
            strcpy(BarretteCourante->NomSignal,Repertoire);
            strcat(BarretteCourante->NomSignal,TempFileName);

            LectureSignal(BarretteCourante,IndiceTempsMax);
        }
    }
    else
    {
        sscanf(LigneCourante,"%i %s",&BarretteCourante->TypeSignal,BarretteCourante->NomSignal);
        LectureSignal(BarretteCourante,IndiceTempsMax);
    }
    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%c",&BarretteCourante->Normale);
    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%i %i",&BarretteCourante->XBar,&BarretteCourante->YBar);
    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%i %i %i %i %i",&BarretteCourante->NbElements,&BarretteCourante->Pitch,
           &BarretteCourante->Largeur,&BarretteCourante->Apodisation,
           &BarretteCourante->Focale);
    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%lf %lf",&BarretteCourante->Direction,
           &BarretteCourante->VitesseMilieuEmission);
}

void LectureParametresBarretteReception(FILE *Fich,char *LigneCourante,struct BARRETTE_RECEPTION *BarretteCourante)
{
    fgets(LigneCourante,500,Fich);
    /*sscanf(LigneCourante,"%s",&BarretteCourante->Nom);bug didier*/
    sscanf(LigneCourante,"%s",BarretteCourante->Nom);
    fgets(LigneCourante,500,Fich);
    //sscanf(LigneCourante,"%c",&BarretteCourante->Normale);

    sscanf(LigneCourante,"%c %i",&BarretteCourante->Normale,&BarretteCourante->ActiveMedium);


//   printf("TOTO \n \n \n");
//   printf("%c \n",BarretteCourante->Normale);
//   printf("%i \n",BarretteCourante->ActiveMedium);

    // exit(1);
//ActiveMedium

    /*    switch(BarretteCourante->Normale)
      {
      case 'C':
      case 'D':
      sscanf(LigneCourante,"%i",&BarretteCourante->ActiveMedium);
      printf("%i \n",BarretteCourante->ActiveMedium);

          break;
      default:
          break;
      }

    */
    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%i %i",
           &BarretteCourante->XBar,&BarretteCourante->YBar);
    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%i %i %i",
           &BarretteCourante->NbElements,&BarretteCourante->Pitch,&BarretteCourante->Largeur);
}
void LectureParametresBarretteSourceFile(FILE *Fich, char *LigneCourante, struct BARRETTE_RECEPTION *BarretteCourante)
{

    FILE *ArraySourceFile=NULL;
    double unused;

    fgets(LigneCourante,500,Fich);
    sscanf(LigneCourante,"%s",BarretteCourante->Nom);

    if((ArraySourceFile=fopen(BarretteCourante->Nom,"rb"))==NULL)
    {
        printf("ArraySourceFile NOT FOUND...\n");
        exit(1);
    }
    else
    {

        fread(&BarretteCourante->Normale,sizeof(char),1,ArraySourceFile);
        fread(&BarretteCourante->NbElements,sizeof(int),1,ArraySourceFile);
        fread(&BarretteCourante->XBar,sizeof(int),1,ArraySourceFile);
        fread(&BarretteCourante->YBar,sizeof(int),1,ArraySourceFile);
        fread(&BarretteCourante->Pitch,sizeof(int),1,ArraySourceFile);
        fread(&BarretteCourante->Largeur,sizeof(int),1,ArraySourceFile);
        fread(&unused,sizeof(double),1,ArraySourceFile);
        fread(&BarretteCourante->DurationPt,sizeof(int),1,ArraySourceFile);

        BarretteCourante->Signaux=(double *) calloc (BarretteCourante->DurationPt*BarretteCourante->NbElements,sizeof(double));
        fread(&unused,sizeof(double),1,ArraySourceFile);
        fread(BarretteCourante->Signaux,sizeof(double),BarretteCourante->DurationPt*BarretteCourante->NbElements,ArraySourceFile);

        fclose(ArraySourceFile );
    }

}
int ValiditeBarretteEmission(struct BARRETTE_EMISSION *BarretteCourante,int LimX,int LimY)
{
    int Test=0;
    switch(BarretteCourante->Normale)
    {
    case 'X':
    case '1':
        Test=((BarretteCourante->XBar < LimX)&&
              ((BarretteCourante->YBar + (BarretteCourante->NbElements-1)*BarretteCourante->Pitch
                + (BarretteCourante->Largeur-1)) < LimY));
        break;
    case 'Y':
    case '2':
        Test=((BarretteCourante->YBar < LimY)&&
              ((BarretteCourante->XBar + (BarretteCourante->NbElements-1)*BarretteCourante->Pitch
                + (BarretteCourante->Largeur-1)) < LimX));
        break;
    default:
        break;
    }
    return Test;
}
int ValiditeBarretteReception(struct BARRETTE_RECEPTION *BarretteCourante,int LimX,int LimY)
{
    int Test=0;
    switch(BarretteCourante->Normale)
    {
    case 'X':
    case '1':
        /*	Default : bug Didier. */
        Test=((BarretteCourante->XBar < LimX) &&
              ((BarretteCourante->YBar + (BarretteCourante->NbElements-1)*BarretteCourante->Pitch
                + (BarretteCourante->Largeur-1)) < LimY));
        break;
    case 'Y':
    case '2':
        Test=((BarretteCourante->YBar < LimY)&&
              ((BarretteCourante->XBar + (BarretteCourante->NbElements-1)*BarretteCourante->Pitch
                + (BarretteCourante->Largeur-1)) < LimX));
        break;
    case 'A':
    case 'B':
    case 'C':
    case 'D':
        Test=1;
    default:
        break;
    }
    return Test;
}

void ActivationBarretteSourceFile(PREC_GRILLE *GrilleCourante, struct BARRETTE_RECEPTION *BarretteCourante,int DimX,int DimY)
{

    int NbElements  =BarretteCourante->NbElements;
    char Normale    =BarretteCourante->Normale;
    int XBar        =BarretteCourante->XBar;
    int YBar        =BarretteCourante->YBar;
    int Pitch       =BarretteCourante->Pitch;
    int Largeur     =BarretteCourante->Largeur;
    double* Signaux =BarretteCourante->Signaux;
    int DurationPt  =BarretteCourante->DurationPt;
    int N,i;
    /* DimX useless for the moment. */
    N=DimX;

    if (IndiceTemps<DurationPt)
    {
        for(N=0; N < NbElements; N++)
        {
            for(i=0; i <Largeur; i++)
            {
                switch(Normale)
                {
                case 'X':
                case '1':
                    if(SourceType==1)
                        GrilleCourante[(XBar)*DimY+(YBar+N*Pitch+i)] +=Signaux[N*DurationPt+IndiceTemps];
                    else if(SourceType==2)
                        GrilleCourante[(XBar)*DimY+(YBar+N*Pitch+i)]=Signaux[N*DurationPt+IndiceTemps];
                    break;
                case 'Y':
                case '2':
                    if(SourceType==1)
                        GrilleCourante[(XBar+N*Pitch+i)*DimY+(YBar)] +=Signaux[N*DurationPt+IndiceTemps];
                    else if(SourceType==2)
                        GrilleCourante[(XBar+N*Pitch+i)*DimY+(YBar)]=Signaux[N*DurationPt+IndiceTemps];
                    break;
                default:
                    break;
                }
            }
        }
    }
    /* Faire ça, puis traiter tous les cas de variables : Txx, Tyy etc... */
}

void ActivationBarretteEmission(PREC_GRILLE *GrilleCourante, struct BARRETTE_EMISSION *BarretteCourante,int DimX, int DimY)
{
    /* DimX useless for the moment */
    /* Didier */
    int N=DimX,i,Retard=0;
    double RetardMaxDeviation=0.0;
    double* CoeffApodisation;
    /* Simplement pour faciliter la relecture du code ... */
    int NbElements  =BarretteCourante->NbElements;
    char Normale    =BarretteCourante->Normale;
    int XBar        =BarretteCourante->XBar;
    int YBar        =BarretteCourante->YBar;
    int Pitch       =BarretteCourante->Pitch;
    int Largeur     =BarretteCourante->Largeur;
    int Apodisation =BarretteCourante->Apodisation;
    int Focale		=BarretteCourante->Focale;
    double Direction=BarretteCourante->Direction;
    double VMilieu  =BarretteCourante->VitesseMilieuEmission;
    double* Signal  =BarretteCourante->Signal;
    int DurationPt  =BarretteCourante->DurationPt;/* didier */
    switch(Normale)
    {
    case 'Y':
    case '2':
        Direction=-Direction; /*Pour rendre plus cohérente la façon d'orienter les directions...*/
        break;
    default:
        break;
    }
    if(VMilieu)
        RetardMaxDeviation=(double)((NbElements-1)*Pitch*H)/VMilieu*sin(2*3.1415/360.0*Direction);

    CoeffApodisation=(double *) calloc (NbElements,sizeof(double));
    switch(Apodisation)
    {
    case 0:
    default:
        for(N=0; N<NbElements; N++)
            CoeffApodisation[N]=1;
        break;
    case 1:
        /* START BUG DIVIDE BY 0 */
        if(NbElements==1)
            CoeffApodisation[0]=0;
        else
        {
            for(N=0; N<NbElements; N++)
                CoeffApodisation[N]=(1-cos(2*3.141593*((double)(N))/((double)(NbElements-1))));
        }
        /*
        for(N=0;N<NbElements;N++)
            CoeffApodisation[N]=(1-cos(2*3.141593*((double)(N))/((double)(NbElements-1))));
        */
        /* END BUG DIVIDE BY 0 */
        break;
    case 2 :
        /* START BUG DIVIDE BY 0 */
        if(NbElements==1)
            CoeffApodisation[0]=1;
        else
        {
            for(N=0; N<NbElements; N++)
                CoeffApodisation[N]=(cos(0.5*3.141593*((double)(N))/((double)(NbElements-1))));
        }
        /*
        for(N=0;N<NbElements;N++)
            CoeffApodisation[N]=(cos(0.5*3.141593*((double)(N))/((double)(NbElements-1))));
        */
        /* END BUG DIVIDE BY 0 */
        break;
    }
    for(N=0; N<NbElements; N++)
    {
        for(i=0; i < Largeur; i++)
        {
            /*Ajout du retard de focalisation (loi exacte, à l'échantillonnage temporel prêt)*/
            if(Focale > 0)
            {
                Retard=(int)((1./DeltaT)*(H/VMilieu*(
                                              sqrt(Focale*Focale+((-(NbElements-1.)/2.+N)*Pitch)*((-(NbElements-1.)/2.+N)*Pitch))
                                              -
                                              sqrt(Focale*Focale+((NbElements-1)/2.*Pitch)*((NbElements-1)/2.*Pitch))
                                          )));
            }
            else
            {
                if(Focale <0)
                {

                    Retard=(int)((-1./DeltaT)*(H/VMilieu*(
                                                   sqrt(Focale*Focale+((-(NbElements-1.)/2.+N)*Pitch)*((-(NbElements-1.)/2.+N)*Pitch))
                                                   -
                                                   sqrt(Focale*Focale)
                                               )));
                }
                else
                    Retard=0;
            }


            /*Ajout du retard de direction (linéaire...)*/
            if(Direction < 0)
            {
                /* START BUG DIVIDE BY 0 */
                if(NbElements==1)
                    Retard +=(int)(1./DeltaT*RetardMaxDeviation);
                else
                    Retard +=(int)(1./DeltaT*RetardMaxDeviation*(1-(double)(N)/(double)(NbElements-1)));
                /*
                Retard +=(int)(1./DeltaT*RetardMaxDeviation*(1-(double)(N)/(double)(NbElements-1)));
                */
                /* END BUG DIVIDE BY 0 */
            }
            else
            {
                /* START BUG DIVIDE BY 0 */
                if(NbElements!=1)
                    Retard +=-(int)(1./DeltaT*RetardMaxDeviation*((double)(N)/(double)(NbElements-1)));
                /*
                Retard +=-(int)(1./DeltaT*RetardMaxDeviation*((double)(N)/(double)(NbElements-1)));
                */
                /* END BUG DIVIDE BY 0 */
            }
            if(((IndiceTemps+Retard)>=0)&&((IndiceTemps+Retard)<IndiceTempsMax))
            {
                switch(Normale)
                {
                case 'X' :
                case '1' :
                    if(SourceType==1)/*20050930*/
                        GrilleCourante[(XBar)*DimY+(YBar+N*Pitch+i)] +=Signal[IndiceTemps+Retard]*CoeffApodisation[N];
                    else if(SourceType==2)
                        if(IndiceTemps < DurationPt)/*20050930*/
                            GrilleCourante[(XBar)*DimY+(YBar+N*Pitch+i)]=Signal[IndiceTemps+Retard]*CoeffApodisation[N];
                    break;
                case 'Y' :
                case '2' :
                    if(SourceType==1)
                        GrilleCourante[(XBar+N*Pitch+i)*DimY+(YBar)] +=Signal[IndiceTemps+Retard]*CoeffApodisation[N];
                    else if(SourceType==2)
                        if(IndiceTemps < DurationPt)/*20050930*/
                            GrilleCourante[(XBar+N*Pitch+i)*DimY+(YBar)]=Signal[IndiceTemps+Retard]*CoeffApodisation[N];
                    break;
                default:
                    break;
                }
            }
        }
    }
    free(CoeffApodisation);
}


/************************ Lectures d'entrées **************************************/

/*lecture du fichier de parametres*/
void LectureParametres(void)
{

    struct BARRETTE_EMISSION *BarretteCouranteEmission;
    struct BARRETTE_RECEPTION *BarretteCouranteReception;

    char ParameterFileName[500];
    char CurrentLine[500];
    char NomFichierSignalPhotoAc[200];

    FILE *Fich=NULL;
    FILE *FichSignalPhotoAc;

    int i=0, IndiceMateriau, Test, NbPtsAlloues, DurationPhotoAcPt;
    double temp=0;

    double CoeffPhotoAc,densite,C11,C12,C22,C33,X11,X12,X22,X33,QModel,VitesseQModel;

    /*Default parameters, corresponding to optional parameters*/
    H=0.1;
    Vmax=1.5;
    Duree=0.0;
    CoeffDeltaT=0.99;
    TypeAbsorption=0;
    IsSkipFileCreated=0;
    IsDoneFileCreated=0;
    ConditionBordXm=0,ConditionBordXp=0,ConditionBordYm=0,ConditionBordYp=0;
    CoeffPML=80,W=40,VmaxPML=1.5;
    NumberSelectedSnapshots=0;


    RcvDwnSmplng=1;
    SourceType=1;
    SnapshotMaxV=0,SnapshotMaxVx=0,SnapshotMaxVy=0,SnapshotMaxTxx=0,SnapshotMaxTyy=0,SnapshotMaxTxy=0, SnapshotMaxEc=0;
    SnapshotV=0,SnapshotVx=0,SnapshotVy=0,SnapshotTxx=0,SnapshotTyy=0,SnapshotTxy=0, SnapshotEc=0;
    SnapshotUx=0;
    SnapshotVx=0;
    PeriodImage=1;
    PeriodRecepteurs=100;
    IndiceTempsAffichageTempsRestant=1;
    NbBarrettesEmissionTxx=0,NbBarrettesEmissionTyy=0,NbBarrettesEmissionTxy=0,NbBarrettesEmissionVx=0,NbBarrettesEmissionVy=0;
    NbBarrettesReceptionTxx=0,NbBarrettesReceptionTyy=0,NbBarrettesReceptionTxy=0,NbBarrettesReceptionVx=0,NbBarrettesReceptionVy=0;
    NbBarrettesSourceFilesTxx=0,NbBarrettesSourceFilesTyy=0,NbBarrettesSourceFilesTxy=0,NbBarrettesSourceFilesVx=0,NbBarrettesSourceFilesVy=0;

    strcpy(ParameterFileName,Repertoire);
    strcat(ParameterFileName,"""Parametres.ini2D");

    if((Fich=fopen(ParameterFileName,"rb"))==NULL)
    {
        strcpy(ParameterFileName,Repertoire);
        strcat(ParameterFileName,"""Parameters.ini2D");

        if((Fich=fopen(ParameterFileName,"rb"))==NULL)
        {

            printf("Parametres.ini2D NOT FOUND ...\n");
            exit(1);
        }
    }

    /*****************************************************************/



    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Vmax    ") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&Vmax);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Grid Step") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&H);
        }
    }



    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Absorption Type") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&TypeAbsorption);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Simulation Length") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&Duree);
        }
    }



    /* for compatibility with previous versions. Fields have been changed...*/
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"DeltaT Coefficient") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&CoeffDeltaT);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"CFL Coefficient") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&CoeffDeltaT);
        }
    }





    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Snapshots Record Period") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&PeriodImage);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Image Record Period") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&PeriodImage);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Data Record Frequency") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&PeriodImage);
            PeriodImage=1./PeriodImage;
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Signals Refresh Period") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&PeriodRecepteurs);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Receivers Record Frequency") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&PeriodRecepteurs);
            PeriodRecepteurs=1./PeriodRecepteurs;
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Time information step period") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&IndiceTempsAffichageTempsRestant);
        }
    }

    /*Parametres déduits des précedents*/
    DeltaT=CoeffDeltaT*H/((sqrt(2.))*Vmax);
    DeltaTsurH=DeltaT/H;
    IndiceTempsMax=(int) (Duree/DeltaT+1);
    IndiceTempsImage=(int)(PeriodImage/DeltaT+1);
    IndiceTempsRecepteurs=(int)(PeriodRecepteurs/DeltaT+1);



    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Number of selected Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NumberSelectedSnapshots);

            SnapshotsTimeSelection=(int *) calloc (NumberSelectedSnapshots,sizeof(int));

            for(i=0; i<NumberSelectedSnapshots; i++)
            {
                fscanf (Fich, "%lf", &temp);
                SnapshotsTimeSelection[i]=(int)floor(temp/DeltaT+0.5);
                printf ("%i \n",SnapshotsTimeSelection[i]);
            }
        }
    }




    /* for compatibility with previous versions. Fields have been changed...*/

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Xm Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordXm);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"X1_low Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordXm);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Xp Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordXp);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"X1_high Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordXp);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Ym Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordYm);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"X2_low Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordYm);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Yp Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordYp);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"X2_high Boundary") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&ConditionBordYp);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Attenuation coeff.in PML") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&CoeffPML);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"PML Efficiency") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&CoeffPML);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"PML Width W") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&W);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"PML Thickness") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&W);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Vmax in PML") !=NULL)
        {
            sscanf(CurrentLine+31,"%lf",&VmaxPML);
        }
    }

    if(TypeAbsorption!=0)
    {

        fseek(Fich,0L,SEEK_SET);
        while(fgets(CurrentLine,500,Fich) !=NULL)
        {
            if(strstr(CurrentLine,"Starts Materials List") !=NULL)
            {
                while(strstr(fgets(CurrentLine,500,Fich),"Ends Materials List")==NULL)
                {
                    sscanf(CurrentLine,"%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                           &IndiceMateriau,&densite,&C11,&C22,&C12,&C33,&X11,&X22,&X12,&X33,&QModel,&VitesseQModel);
                    TableLegerete[(PREC_CARTE)IndiceMateriau]=1./densite;
                    TableC11[(PREC_CARTE)IndiceMateriau]=C11;
                    TableC22[(PREC_CARTE)IndiceMateriau]=C22;
                    TableC12[(PREC_CARTE)IndiceMateriau]=C12;
                    TableC33[(PREC_CARTE)IndiceMateriau]=C33;
                    TableX11[(PREC_CARTE)IndiceMateriau]=X11;
                    TableX22[(PREC_CARTE)IndiceMateriau]=X22;
                    TableX12[(PREC_CARTE)IndiceMateriau]=X12;
                    TableX33[(PREC_CARTE)IndiceMateriau]=X33;
                    TableQModel[(PREC_CARTE)IndiceMateriau]=exp(-log(10)/20*QModel*VitesseQModel*DeltaT);
                }

            }
        }
    }

    if(TypeAbsorption==0)
    {
        fseek(Fich,0L,SEEK_SET);
        while(fgets(CurrentLine,500,Fich) !=NULL)
        {
            if(strstr(CurrentLine,"Starts Materials List") !=NULL)
            {
                while(strstr(fgets(CurrentLine,500,Fich),"Ends Materials List")==NULL)
                {
                    sscanf(CurrentLine,"%i %lf %lf %lf %lf %lf",
                           &IndiceMateriau,&densite,&C11,&C22,&C12,&C33);
                    TableLegerete[(PREC_CARTE)IndiceMateriau]=1./densite;
                    TableC11[(PREC_CARTE)IndiceMateriau]=C11;
                    TableC22[(PREC_CARTE)IndiceMateriau]=C22;
                    TableC12[(PREC_CARTE)IndiceMateriau]=C12;
                    TableC33[(PREC_CARTE)IndiceMateriau]=C33;
                    TableX11[(PREC_CARTE)IndiceMateriau]=0.;
                    TableX22[(PREC_CARTE)IndiceMateriau]=0.;
                    TableX12[(PREC_CARTE)IndiceMateriau]=0.;
                    TableX33[(PREC_CARTE)IndiceMateriau]=0.;
                    TableQModel[(PREC_CARTE)IndiceMateriau]=1.;
                }
            }
        }
    }

    /*****************************************************************/
    /********************* Parametres optionnels *********************/
    /*****************************************************************/

    /* Options for creating skip/done file*/

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Is skip file created") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&IsSkipFileCreated);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Is done file created") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&IsDoneFileCreated);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Time Reversal mode") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&Is_TR_mode_with_initial_condition);
        }
    }



    /* Options for recording snapshots */
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotV);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotV);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V Data") !=NULL)  /* for compatibility*/
        {
            sscanf(CurrentLine+31,"%i",&SnapshotV);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Ekinetic Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotEc);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Ec Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotEc);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Ec Data") !=NULL) /* for compatibility*/
        {
            sscanf(CurrentLine+31,"%i",&SnapshotEc);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Vx Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotVx);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V1 Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotVx);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Vx Data") !=NULL) /* for compatibility*/
        {
            sscanf(CurrentLine+31,"%i",&SnapshotVx);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Vy Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotVy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V2 Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotVy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Vy Data") !=NULL) /* for compatibility*/
        {
            sscanf(CurrentLine+31,"%i",&SnapshotVy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Ux Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotUx);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record U1 Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotUx);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Uy Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotUy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record U2 Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotUy);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Txx Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTxx);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T11 Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTxx);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Txx Data") !=NULL) /* for compatibility*/
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTxx);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Tyy Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTyy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T22 Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTyy);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Tyy Data") !=NULL) /* for compatibility*/
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTyy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Txy Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTxy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T12 Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTxy);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Txy Data") !=NULL) /* for compatibility*/
        {
            sscanf(CurrentLine+31,"%i",&SnapshotTxy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Curl Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotCurl);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Curl Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotCurl);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Div Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotDiv);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Div Images") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotDiv);
        }
    }

    /* Options for recording maximum images */
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxV);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Ekinetic_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxEc);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V1_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxVx);
        }
    }




    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V2_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxVy);
        }
    }



    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T11_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTxx);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T22_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTyy);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T12_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTxy);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Curl_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxCurl);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Div_max Snapshots") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxDiv);
        }
    }




    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxV);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Ec MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxEc);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Vx MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxVx);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V1 MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxVx);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Vy MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxVy);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record V2 MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxVy);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Txx MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTxx);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Tyy MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTyy);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Txy MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTxy);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T11 MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTxx);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T22 MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTyy);
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record T12 MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxTxy);
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Curl MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxCurl);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Record Div MaxImage") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SnapshotMaxDiv);
        }
    }




    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Receivers Downsampling") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&RcvDwnSmplng);
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Images Folder") !=NULL)
        {
            sscanf(CurrentLine+31,"%s",ImageFolder);
            FoundImageFolder=1;
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Type of Source") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&SourceType);
        }
    }





    /*****************************************************************/
    /**********Lecture des parametres de barrette d'EMISSION *********/
    /*****************************************************************/


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Txx Emitting Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionTxx);
            if(NbBarrettesEmissionTxx)
                EmissionTxx=(struct BARRETTE_EMISSION *) calloc(NbBarrettesEmissionTxx,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionTxx; i++)
            {
                BarretteCouranteEmission=&EmissionTxx[i];
                LectureParametresBarretteEmission(Fich,CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X,Y);
                if(!Test)
                {
                    printf("Pb on emitters array T11 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T11 Emitters Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionTxx);
            if(NbBarrettesEmissionTxx)
                EmissionTxx=(struct BARRETTE_EMISSION *) calloc(NbBarrettesEmissionTxx,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionTxx; i++)
            {
                BarretteCouranteEmission=&EmissionTxx[i];
                LectureParametresBarretteEmission(Fich,CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X,Y);
                if(!Test)
                {
                    printf("Pb on emitters array T11 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Tyy Emitting Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionTyy);
            if(NbBarrettesEmissionTyy)
                EmissionTyy=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionTyy,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionTyy; i++)
            {
                BarretteCouranteEmission=&EmissionTyy[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X,Y);
                if (!Test)
                {
                    printf("Pb on emitters array T22 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T22 Emitters Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionTyy);
            if(NbBarrettesEmissionTyy)
                EmissionTyy=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionTyy,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionTyy; i++)
            {
                BarretteCouranteEmission=&EmissionTyy[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X,Y);
                if (!Test)
                {
                    printf("Pb on emitters array T22 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }



    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Txy Emitting Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionTxy);
            if(NbBarrettesEmissionTxy)
                EmissionTxy=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionTxy,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionTxy; i++)
            {
                BarretteCouranteEmission=&EmissionTxy[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X+1,Y+1);
                if(!Test)
                {
                    printf("Pb on emitters array T12 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T12 Emitters Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionTxy);
            if(NbBarrettesEmissionTxy)
                EmissionTxy=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionTxy,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionTxy; i++)
            {
                BarretteCouranteEmission=&EmissionTxy[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X+1,Y+1);
                if(!Test)
                {
                    printf("Pb on emitters array T12 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }




    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Vx Emitting Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionVx);
            if(NbBarrettesEmissionVx)
                EmissionVx=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionVx,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionVx; i++)
            {
                BarretteCouranteEmission=&EmissionVx[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X+1,Y);
                if(!Test)
                {
                    printf("Pb on emitters array V1 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of V1 Emitters Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionVx);
            if(NbBarrettesEmissionVx)
                EmissionVx=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionVx,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionVx; i++)
            {
                BarretteCouranteEmission=&EmissionVx[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X+1,Y);
                if(!Test)
                {
                    printf("Pb on emitters array V1 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Vy Emitting Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionVy);
            if(NbBarrettesEmissionVy)
                EmissionVy=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionVy,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionVy; i++)
            {
                BarretteCouranteEmission=&EmissionVy[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X,Y+1);
                if(!Test)
                {
                    printf("Pb on emitters array V2 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of V2 Emitters Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesEmissionVy);
            if(NbBarrettesEmissionVy)
                EmissionVy=(struct BARRETTE_EMISSION*) calloc(NbBarrettesEmissionVy,sizeof(struct BARRETTE_EMISSION));
            for(i=0; i < NbBarrettesEmissionVy; i++)
            {
                BarretteCouranteEmission=&EmissionVy[i];
                LectureParametresBarretteEmission(Fich, CurrentLine,BarretteCouranteEmission);
                Test=ValiditeBarretteEmission(BarretteCouranteEmission,X,Y+1);
                if(!Test)
                {
                    printf("Pb on emitters array V2 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }



    /*****************************************************************/
    /**********Lecture des fichiers source de barrette d'EMISSION *********/
    /*****************************************************************/

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T11 Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesTxx);
            if(NbBarrettesSourceFilesTxx)
                BarrettesSourceFilesTxx=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesTxx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesTxx; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesTxx[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on emission source file T11 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T22 Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesTyy);
            if(NbBarrettesSourceFilesTyy)
                BarrettesSourceFilesTyy=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesTyy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesTyy; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesTyy[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on emission source file T22 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T12 Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesTxy);
            if(NbBarrettesSourceFilesTxy)
                BarrettesSourceFilesTxy=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesTxy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesTxy; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesTxy[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y+1);
                if(!Test)
                {
                    printf("Pb on emission source file T12 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of V1 Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesVx);
            if(NbBarrettesSourceFilesVx)
                BarrettesSourceFilesVx=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesVx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesVx; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesVx[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y);
                if(!Test)
                {
                    printf("Pb on emission source file V1 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of V2 Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesVy);
            if(NbBarrettesSourceFilesVy)
                BarrettesSourceFilesVy=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesVy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesVy; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesVy[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y+1);
                if(!Test)
                {
                    printf("Pb on emission source file V2 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Txx Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesTxx);
            if(NbBarrettesSourceFilesTxx)
                BarrettesSourceFilesTxx=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesTxx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesTxx; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesTxx[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on emission source file T11 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Tyy Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesTyy);
            if(NbBarrettesSourceFilesTyy)
                BarrettesSourceFilesTyy=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesTyy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesTyy; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesTyy[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on emission source file T22 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Txy Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesTxy);
            if(NbBarrettesSourceFilesTxy)
                BarrettesSourceFilesTxy=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesTxy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesTxy; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesTxy[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y+1);
                if(!Test)
                {
                    printf("Pb on emission source file T12 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Vx Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesVx);
            if(NbBarrettesSourceFilesVx)
                BarrettesSourceFilesVx=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesVx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesVx; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesVx[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y);
                if(!Test)
                {
                    printf("Pb on emission source file V1 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Vy Array Source Files") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesSourceFilesVy);
            if(NbBarrettesSourceFilesVy)
                BarrettesSourceFilesVy=(struct BARRETTE_RECEPTION *) calloc(NbBarrettesSourceFilesVy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesSourceFilesVy; i++)
            {
                BarretteCouranteReception=&BarrettesSourceFilesVy[i];
                LectureParametresBarretteSourceFile(Fich,CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y+1);
                if(!Test)
                {
                    printf("Pb on emission source file V2 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
            }
        }
    }


    /*****************************************************************/
    /**********Lecture des parametres de barrette de RECEPTION *******/
    /*****************************************************************/

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Txx Receiving Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionTxx);
            if(NbBarrettesReceptionTxx)
                ReceptionTxx=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionTxx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionTxx; i++)
            {
                BarretteCouranteReception=&ReceptionTxx[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on receivers array T11 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }

                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T11 Receivers Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionTxx);

            if(NbBarrettesReceptionTxx)

                ReceptionTxx=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionTxx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionTxx; i++)
            {
                BarretteCouranteReception=&ReceptionTxx[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on receivers array T11 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }

                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }





    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Tyy Receiving Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionTyy);
            if(NbBarrettesReceptionTyy)
                ReceptionTyy=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionTyy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionTyy; i++)
            {
                BarretteCouranteReception=&ReceptionTyy[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on receivers array T22 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T22 Receivers Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionTyy);
            if(NbBarrettesReceptionTyy)
                ReceptionTyy=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionTyy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionTyy; i++)
            {
                BarretteCouranteReception=&ReceptionTyy[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y);
                if(!Test)
                {
                    printf("Pb on receivers array T22 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }




    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Txy Receiving Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionTxy);
            if(NbBarrettesReceptionTxy)
                ReceptionTxy=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionTxy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionTxy; i++)
            {
                BarretteCouranteReception=&ReceptionTxy[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y+1);
                if(!Test)
                {
                    printf("Pb on receivers array T12 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }
    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of T12 Receivers Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionTxy);
            if(NbBarrettesReceptionTxy)
                ReceptionTxy=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionTxy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionTxy; i++)
            {
                BarretteCouranteReception=&ReceptionTxy[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y+1);
                if(!Test)
                {
                    printf("Pb on receivers array T12 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }





    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Vx Receiving Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionVx);
            if(NbBarrettesReceptionVx)
                ReceptionVx=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionVx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionVx; i++)
            {
                BarretteCouranteReception=&ReceptionVx[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y);
                if(!Test)
                {
                    printf("Pb on receivers array V1 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }


    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of V1 Receivers Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionVx);
            if(NbBarrettesReceptionVx)
                ReceptionVx=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionVx,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionVx; i++)
            {
                BarretteCouranteReception=&ReceptionVx[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X+1,Y);
                if(!Test)
                {
                    printf("Pb on receivers array V1 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }




    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of Vy Receiving Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionVy);
            if(NbBarrettesReceptionVy)
                ReceptionVy=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionVy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionVy; i++)
            {
                BarretteCouranteReception=&ReceptionVy[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y+1);
                if(!Test)
                {
                    printf("Pb on receivers array V2 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }

    fseek(Fich,0L,SEEK_SET);
    while(fgets(CurrentLine,500,Fich) !=NULL)
    {
        if(strstr(CurrentLine,"Nb of V2 Receivers Array") !=NULL)
        {
            sscanf(CurrentLine+31,"%i",&NbBarrettesReceptionVy);
            if(NbBarrettesReceptionVy)
                ReceptionVy=(struct BARRETTE_RECEPTION*) calloc(NbBarrettesReceptionVy,sizeof(struct BARRETTE_RECEPTION));
            for(i=0; i < NbBarrettesReceptionVy; i++)
            {
                BarretteCouranteReception=&ReceptionVy[i];
                LectureParametresBarretteReception(Fich, CurrentLine,BarretteCouranteReception);
                Test=ValiditeBarretteReception(BarretteCouranteReception,X,Y+1);
                if(!Test)
                {
                    printf("Pb on receivers array V2 number %i\n",i+1);
                    LiberationMemoire();
                    exit(1);
                }
                NbPtsAlloues=(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng))*(BarretteCouranteReception->NbElements);
                BarretteCouranteReception->Signaux=(double *) calloc(NbPtsAlloues,sizeof(double));
            }
        }
    }
    /* ModePhotoAcPhotoAc  is turned to 1 if a compatible photoacoustic
    map is encountered*/
    if(ModePhotoAc)
    {
        fseek(Fich, 0L, SEEK_SET);
        while(strstr(fgets(CurrentLine,500,Fich),"Starts Coefficient List")==NULL);
        while(strstr(fgets(CurrentLine,500,Fich),"Ends Coefficient List")==NULL)
        {
            sscanf(CurrentLine,"%i %lf",&IndiceMateriau,&CoeffPhotoAc);
            TableCoeffPhotoAc[(PREC_CARTE)IndiceMateriau]=CoeffPhotoAc;
        }
        fseek(Fich, 0L, SEEK_SET);
        while(fgets(CurrentLine,500,Fich) !=NULL)
        {
            if(strstr(CurrentLine,"Photoacoustic Signal") !=NULL)
            {
                sscanf(CurrentLine+31,"%s",NomFichierSignalPhotoAc);
            }
        }
        SignalPhotoAc=(double *) calloc (IndiceTempsMax,sizeof(double));
        if((FichSignalPhotoAc=fopen(NomFichierSignalPhotoAc, "rb" ))==NULL)
        {
            printf("Problem opening photoacoustic signal file %s: file not found!\n",NomFichierSignalPhotoAc);
            LiberationMemoire();
            exit(1);
        }
        fread(&DurationPhotoAcPt,sizeof(int),1,FichSignalPhotoAc );
        fread(SignalPhotoAc,sizeof(double),DurationPhotoAcPt,FichSignalPhotoAc);
        fclose(FichSignalPhotoAc);
    }
    fclose(Fich);
}

void LectureConditionsInitiales(char *NomFichierCondIni,PREC_GRILLE *GrilleLue,int DimX, int DimY)
{
    int Entete,i;
    FILE *Fich=NULL;

    Fich=fopen(NomFichierCondIni,"rb");
    Entete=2*sizeof(int)+3*sizeof(double);
    fseek(Fich,Entete, SEEK_SET);
    for(i=0; i < DimX*DimY; i++)
    {
        fread(GrilleLue+i,sizeof(float),1,Fich);
    }
    fclose(Fich);
}

void LectureSignal(struct BARRETTE_EMISSION *BarretteCourante, int NbPtsSignal)
{
    FILE *Fich=NULL;
    if((Fich=fopen(BarretteCourante->NomSignal,"rb"))==NULL)
    {
        printf("Problem opening signal file %s :file not found!\n",BarretteCourante->NomSignal);
        LiberationMemoire();
        exit(1);
    }
    /*Lecture des dimensions de la simulation dans les 2 premiers octets de carte.map*/
    fread(&BarretteCourante->DurationPt,sizeof(int),1,Fich);
    if(NbPtsSignal < BarretteCourante->DurationPt)
    {
        printf("Emission signals will be truncated, as longer than simulation time...\n");
        BarretteCourante->DurationPt=NbPtsSignal;
    }
    /*Lecture et fermeture fichier du signal*/
    fread(BarretteCourante->Signal,sizeof(double),BarretteCourante->DurationPt,Fich);
    fclose(Fich);
}

/************************ Ecriture de fichiers *************************************/

void SauvegardeRecepteurs(void)
{
    int NumeroBarrette;
    FILE *Fich=NULL;
    char NomFich[200];
    char NomComplet[200];
    int NbPoints;
    double DeltaT_Rcv;
    int IndiceTempsMaxRcv;
    DeltaT_Rcv=DeltaT*RcvDwnSmplng;
    IndiceTempsMaxRcv=1 + ((int)((IndiceTempsMax - 1)/RcvDwnSmplng));
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionTxx; NumeroBarrette++)
    {
        sprintf(NomComplet,Repertoire);
        /*sprintf(NomFich,"");*/ /*didier 5x*/
        /*sprintf(NomFich,"%s",""); aurait été correct*/

        strcpy(NomFich,"");
        strcat(NomFich,ReceptionTxx[NumeroBarrette].Nom);
        strcat(NomFich,"_T11.rcv2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet, "w+b" );
        fwrite(&ReceptionTxx[NumeroBarrette].Normale,sizeof(char),1,Fich);
        fwrite(&ReceptionTxx[NumeroBarrette].NbElements,sizeof(int),1,Fich);
        fwrite(&ReceptionTxx[NumeroBarrette].XBar,sizeof(int),1,Fich);
        fwrite(&ReceptionTxx[NumeroBarrette].YBar,sizeof(int),1,Fich);
        fwrite(&ReceptionTxx[NumeroBarrette].Pitch,sizeof(int),1,Fich);
        fwrite(&ReceptionTxx[NumeroBarrette].Largeur,sizeof(int),1,Fich);
        fwrite(&H,sizeof(double),1,Fich);
        fwrite(&IndiceTempsMaxRcv,sizeof(int),1,Fich);
        fwrite(&DeltaT_Rcv,sizeof(double),1,Fich);
        NbPoints=ReceptionTxx[NumeroBarrette].NbElements*(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng));
        fwrite(ReceptionTxx[NumeroBarrette].Signaux,sizeof(double),NbPoints,Fich);
        fclose(Fich);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionTyy; NumeroBarrette++)
    {
        sprintf(NomComplet,Repertoire);
        strcpy(NomFich,"");
        strcat(NomFich,ReceptionTyy[NumeroBarrette].Nom);
        strcat(NomFich,"_T22.rcv2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        fwrite(&ReceptionTyy[NumeroBarrette].Normale,sizeof(char),1,Fich);
        fwrite(&ReceptionTyy[NumeroBarrette].NbElements,sizeof(int),1,Fich);
        fwrite(&ReceptionTyy[NumeroBarrette].XBar,sizeof(int),1,Fich);
        fwrite(&ReceptionTyy[NumeroBarrette].YBar,sizeof(int),1,Fich);
        fwrite(&ReceptionTyy[NumeroBarrette].Pitch,sizeof(int),1,Fich);
        fwrite(&ReceptionTyy[NumeroBarrette].Largeur,sizeof(int),1,Fich);
        fwrite(&H,sizeof(double),1,Fich);
        fwrite(&IndiceTempsMaxRcv,sizeof(int),1,Fich);
        fwrite(&DeltaT_Rcv,sizeof(double),1,Fich);
        NbPoints=ReceptionTyy[NumeroBarrette].NbElements*(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng));
        fwrite(ReceptionTyy[NumeroBarrette].Signaux,sizeof(double),NbPoints,Fich);
        fclose(Fich);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionTxy; NumeroBarrette++)
    {
        sprintf(NomComplet,Repertoire);
        strcpy(NomFich,"");
        strcat(NomFich,ReceptionTxy[NumeroBarrette].Nom);
        strcat(NomFich,"_T12.rcv2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        fwrite(&ReceptionTxy[NumeroBarrette].Normale,sizeof(char),1,Fich);
        fwrite(&ReceptionTxy[NumeroBarrette].NbElements,sizeof(int),1,Fich);
        fwrite(&ReceptionTxy[NumeroBarrette].XBar,sizeof(int),1,Fich);
        fwrite(&ReceptionTxy[NumeroBarrette].YBar,sizeof(int),1,Fich);
        fwrite(&ReceptionTxy[NumeroBarrette].Pitch,sizeof(int),1,Fich);
        fwrite(&ReceptionTxy[NumeroBarrette].Largeur,sizeof(int),1,Fich);
        fwrite(&H,sizeof(double),1,Fich);
        fwrite(&IndiceTempsMaxRcv,sizeof(int),1,Fich);
        fwrite(&DeltaT_Rcv,sizeof(double),1,Fich);
        NbPoints=ReceptionTxy[NumeroBarrette].NbElements*(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng));
        fwrite(ReceptionTxy[NumeroBarrette].Signaux,sizeof(double),NbPoints,Fich);
        fclose(Fich);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionVx; NumeroBarrette++)
    {
        sprintf(NomComplet,Repertoire);
        strcpy(NomFich,"");
        strcat(NomFich,ReceptionVx[NumeroBarrette].Nom);
        strcat(NomFich,"_V1.rcv2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        fwrite(&ReceptionVx[NumeroBarrette].Normale,sizeof(char),1,Fich);
        fwrite(&ReceptionVx[NumeroBarrette].NbElements,sizeof(int),1,Fich);
        fwrite(&ReceptionVx[NumeroBarrette].XBar,sizeof(int),1,Fich);
        fwrite(&ReceptionVx[NumeroBarrette].YBar,sizeof(int),1,Fich);
        fwrite(&ReceptionVx[NumeroBarrette].Pitch,sizeof(int),1,Fich);
        fwrite(&ReceptionVx[NumeroBarrette].Largeur,sizeof(int),1,Fich);
        fwrite(&H,sizeof(double),1,Fich);
        fwrite(&IndiceTempsMaxRcv,sizeof(int),1,Fich);
        fwrite(&DeltaT_Rcv,sizeof(double),1,Fich);
        NbPoints=ReceptionVx[NumeroBarrette].NbElements*(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng));
        fwrite(ReceptionVx[NumeroBarrette].Signaux,sizeof(double),NbPoints,Fich);
        fclose(Fich);
    }
    for (NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionVy; NumeroBarrette++)
    {
        sprintf(NomComplet,Repertoire);
        strcpy(NomFich,"");
        strcat(NomFich,ReceptionVy[NumeroBarrette].Nom);
        strcat(NomFich,"_V2.rcv2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        fwrite(&ReceptionVy[NumeroBarrette].Normale,sizeof(char),1,Fich);
        fwrite(&ReceptionVy[NumeroBarrette].NbElements,sizeof(int),1,Fich);
        fwrite(&ReceptionVy[NumeroBarrette].XBar,sizeof(int),1,Fich);
        fwrite(&ReceptionVy[NumeroBarrette].YBar,sizeof(int),1,Fich);
        fwrite(&ReceptionVy[NumeroBarrette].Pitch,sizeof(int),1,Fich);
        fwrite(&ReceptionVy[NumeroBarrette].Largeur,sizeof(int),1,Fich);
        fwrite(&H,sizeof(double),1,Fich);
        fwrite(&IndiceTempsMaxRcv,sizeof(int),1,Fich);
        fwrite(&DeltaT_Rcv,sizeof(double),1,Fich);
        NbPoints=ReceptionVy[NumeroBarrette].NbElements*(1+(int)((IndiceTempsMax-1)/RcvDwnSmplng));
        fwrite(ReceptionVy[NumeroBarrette].Signaux,sizeof(double),NbPoints,Fich);
        fclose(Fich);
    }
}
void SauvegardeSnapshot(int NumeroImageCourante)
{
    FILE *Fich=NULL;
    char NomFich[200];
    char NomComplet[200];
    int DimX,DimY;
    int i,j;
    float Temp;
    double when, Legerete;
    double DerivH=1./H;

    /*	 Sauvegarde en float*/
    if(SnapshotTyy)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""T22_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(Tyy[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotTxx)
    {


        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);


        sprintf(NomFich,"""T11_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(Txx[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotTxy)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""T12_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X + 1;
            DimY=Y + 1;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(Txy[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotVx)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""V1_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X + 1;
            DimY=Y;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(Vx[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotUx)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""U1_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X + 1;
            DimY=Y;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(Ux[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotVy)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""V2_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y + 1;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(Vy[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotUy)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""U2_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y + 1;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(Uy[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotV)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""V_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX; i++)
            {
                for(j=0; j < DimY; j++)
                {
                    //       Temp=(float)sqrt(Vy[i*(Y+1)+j]*Vy[i*(Y+1)+j]+Vx[i*Y+j]*Vx[i*Y+j]);
                    Temp=(float)sqrt(0.25*(Vy[i*(Y+1)+j]+Vy[i*(Y+1)+j+1])*(Vy[i*(Y+1)+j]+Vy[i*(Y+1)+j+1])+0.25*(Vx[i*Y+j]+Vx[(i+1)*Y+j])*(Vx[i*Y+j]+Vx[(i+1)*Y+j]));
                    fwrite(&Temp,sizeof(float),1,Fich);
                }
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotEc)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""Ec_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX; i++)
            {
                for(j=0; j < DimY; j++)
                {
                    Legerete=(double)(TableLegerete[Carte[i*Y+j]]);
                    //Temp=(float)(0.5*(Vy[i*(Y+1)+j]*Vy[i*(Y+1)+j]+Vx[i*Y+j]*Vx[i*Y+j])/Legerete);
                    Temp=(float)((0.125*(Vy[i*(Y+1)+j]+Vy[i*(Y+1)+j+1])*(Vy[i*(Y+1)+j]+Vy[i*(Y+1)+j+1])+
                                  0.125*(Vx[i*Y+j]+Vx[(i+1)*Y+j])*(Vx[i*Y+j]+Vx[(i+1)*Y+j]))
                                 /Legerete);
                    fwrite(&Temp,sizeof(float),1,Fich);
                }
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotDiv)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""Div_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX; i++)
            {
                for(j=0; j < DimY; j++)
                {
                    Temp=(float)(DerivH*((Vx[(i+1)*Y+j]-Vx[(i)*Y+j])+(Vy[i*(Y+1)+j+1]-Vy[i*(Y+1)+j])));
                    fwrite(&Temp,sizeof(float),1,Fich);
                }
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotCurl)
    {
        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);

        sprintf(NomFich,"""Curl_%0.3i.snp2D",NumeroImageCourante);
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X - 1;
            DimY=Y - 1;
            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=((double)IndiceTemps)*DeltaT;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX; i++)
            {
                for(j=0; j < DimY; j++)
                {
                    /*Temp=(float)(DerivH*((Vy[(i+1)*(Y+1)+j]-Vy[i*(Y+1)+j])-(Vx[(i)*Y+j+1]-Vx[(i)*Y+j])));*/
                    Temp=(float)(DerivH*((Vy[(i+1)*(Y+1)+j+1]-Vy[i*(Y+1)+j+1])-(Vx[(i+1)*Y+j+1]-Vx[(i+1)*Y+j])));
                    fwrite(&Temp,sizeof(float),1,Fich);
                }
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
}

void RecordMaxSnapshot(void)
{
    FILE *Fich=NULL;
    char NomFich[200];
    char NomComplet[200];
    int DimX,DimY;
    int i;
    float Temp;
    double when;

    if(SnapshotMaxTxx)
    {


        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxT11.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxTxx[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxTyy)
    {


        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxT22.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxTyy[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxTxy)
    {


        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxT12.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X + 1;
            DimY=Y + 1;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxTxy[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxVx)
    {


        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxV1.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X + 1;
            DimY=Y;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxVx[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxVy)
    {

        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxV2.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y + 1;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxVy[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxV)
    {

        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxV.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxV[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxEc)
    {

        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxEc.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxEc[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxCurl)
    {

        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxCurl.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X - 1;
            DimY=Y - 1;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxCurl[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
    if(SnapshotMaxDiv)
    {

        sprintf(NomComplet,Repertoire);
        if (FoundImageFolder==1)
            strcat(NomComplet,ImageFolder);
        sprintf(NomFich,"""MaxDiv.snp2D");
        strcat(NomComplet,NomFich);
        Fich=fopen(NomComplet,"w+b");
        if(Fich!=NULL)
        {
            DimX=X;
            DimY=Y;

            fwrite(&DimX,sizeof(int),1,Fich);
            fwrite(&DimY,sizeof(int),1,Fich);
            when=(double)0;
            fwrite(&when,sizeof(double),1,Fich);
            fwrite(&H,sizeof(double),1,Fich);
            fwrite(&DeltaT,sizeof(double),1,Fich);
            for(i=0; i < DimX*DimY; i++)
            {
                Temp=(float)(MaxDiv[i]);
                fwrite(&Temp,sizeof(float),1,Fich);
            }
            fclose(Fich);
        }
        else
        {
            printf("Problem with directory ImageData\n");
            LiberationMemoire();
            exit(1);
        }
    }
}
void EnregistrementBarretteReception(PREC_GRILLE *GrilleCourante,struct BARRETTE_RECEPTION *BarretteCourante, int DimX,int DimY)
{
    /* DimX useless for the moment. */
    int N=DimX,i;
    int NbElements=BarretteCourante->NbElements;
    int XBar=BarretteCourante->XBar;
    int YBar=BarretteCourante->YBar;
    int Largeur=BarretteCourante->Largeur;
    int Pitch=BarretteCourante->Pitch;
    double *Signaux=BarretteCourante->Signaux;
    int IndiceTempsMaxRcv=1 + ((int)((IndiceTempsMax-1)/RcvDwnSmplng));
    int IndiceTempsRcv=IndiceTemps/RcvDwnSmplng;
    int count=0;
    /* printf("%i",IndiceTempsRcv); */
    switch(BarretteCourante->Normale)
    {
    case 'X' :
    case '1' :

        for(N=0; N < NbElements; N++)
        {
            for(i=0; i < Largeur; i++)
            {
                Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv] +=
                    GrilleCourante[(XBar)*DimY+(YBar+N*Pitch+i)];
            }
        }
        break;
    case 'Y' :
    case '2':
        for(N=0; N < NbElements; N++)
        {
            for(i=0; i < Largeur; i++)
            {
                Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv] +=
                    GrilleCourante[(XBar+N*Pitch+i)*DimY+(YBar)];
            }
        }
        break;
    case 'A' :
        for(N=0; N < NbElements; N++)
        {
            for(i=0; i < Largeur; i++)
            {
                Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv] +=
                    GrilleCourante[(XBar+i)*DimY+(YBar+N*Pitch)];
            }
        }
        break;
    case 'B' :
        for(N=0; N < NbElements; N++)
        {
            for(i=0; i < Largeur; i++)
            {
                Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv] +=
                    GrilleCourante[(XBar+N*Pitch)*DimY+(YBar+i)];
            }
        }
        break;
    case 'C' :
        for(N=0; N < NbElements; N++)
        {
            count=0;
            for(i=0; i < Largeur; i++)

            {
                if (Carte[(XBar+i)*DimY+(YBar+N*Pitch)]==(BarretteCourante->ActiveMedium))
                {
                    count+=1;
                    Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv] +=
                        GrilleCourante[(XBar+i)*DimY+(YBar+N*Pitch)];

                }
            }
            if(count>0)
            {
                Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv]/=(double)count;
                //   printf("%i ",count);
            }
        }
        break;
    case 'D' :
        for(N=0; N < NbElements; N++)
        {
            count=0;
            for(i=0; i < Largeur; i++)

            {
                if (Carte[(XBar+N*Pitch)*DimY+(YBar+i)]==(BarretteCourante->ActiveMedium))
                {
                    count+=1;
                    Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv] +=
                        GrilleCourante[(XBar+N*Pitch)*DimY+(YBar+i)];

                }
            }
            if(count>0)
            {
                Signaux[N*IndiceTempsMaxRcv+IndiceTempsRcv]/=count;
                //   printf("%i ",count);
            }
        }
        break;
    default:
        break;
    }
}



/************************ gestion mémoire ******************************************/

void AllocationGrilles(void)
{
    int DimensionMax;

    Carte=(PREC_CARTE *) calloc(X*Y,sizeof(PREC_CARTE));
    AttenuationPML=(double*) calloc((2*W+1),sizeof(double));

    if (ModePhotoAc)
    {
        CartePhotoAc=(PREC_CARTE *) calloc(X*Y,sizeof(PREC_CARTE));
    }

    /***************************************************************************/
    /*********** Allocation des grilles du bloc principal **********************/
    /***************************************************************************/
    Vx    =(PREC_GRILLE *) calloc((X+1)*Y,sizeof(PREC_GRILLE));
    Vy    =(PREC_GRILLE *) calloc((X)*(Y+1),sizeof(PREC_GRILLE));
    Txx   =(PREC_GRILLE *) calloc(X*Y,sizeof(PREC_GRILLE));
    Tyy   =(PREC_GRILLE *) calloc(X*Y,sizeof(PREC_GRILLE));
    Txy   =(PREC_GRILLE *) calloc((X+1)*(Y+1),sizeof(PREC_GRILLE));

    if(SnapshotUx==1)
        Ux   =(PREC_GRILLE *) calloc((X+1)*Y,sizeof(PREC_GRILLE));

    if(SnapshotUy==1)
        Uy   =(PREC_GRILLE *) calloc((X)*(Y+1),sizeof(PREC_GRILLE));

    if(SnapshotMaxTxx==1)
        MaxTxx   =(PREC_GRILLE *) calloc(X*Y,sizeof(PREC_GRILLE));

    if(SnapshotMaxTyy==1)
        MaxTyy   =(PREC_GRILLE *) calloc(X*Y,sizeof(PREC_GRILLE));

    if(SnapshotMaxTxy==1)
        MaxTxy   =(PREC_GRILLE *) calloc((X+1)*(Y+1),sizeof(PREC_GRILLE));

    if(SnapshotMaxVx==1)
        MaxVx    =(PREC_GRILLE *) calloc((X+1)*Y,sizeof(PREC_GRILLE));

    if(SnapshotMaxVy==1)
        MaxVy    =(PREC_GRILLE *) calloc((X)*(Y+1),sizeof(PREC_GRILLE));

    if(SnapshotMaxV==1)
        MaxV   =(PREC_GRILLE *) calloc(X*Y,sizeof(PREC_GRILLE));

    if(SnapshotMaxEc==1)
        MaxEc   =(PREC_GRILLE *) calloc(X*Y,sizeof(PREC_GRILLE));

    if(SnapshotMaxCurl==1)
        MaxCurl   =(PREC_GRILLE *) calloc((X-1)*(Y-1),sizeof(PREC_GRILLE));

    if(SnapshotMaxDiv==1)
        MaxDiv   =(PREC_GRILLE *) calloc(X*Y,sizeof(PREC_GRILLE));

    if(TypeAbsorption==1)
    {
        Vx_Old=(PREC_GRILLE *) calloc((X+1)*Y,sizeof(PREC_GRILLE));
        Vy_Old=(PREC_GRILLE *) calloc((X)*(Y+1),sizeof(PREC_GRILLE));
    }
    /***************************************************************************/
    /****************** Allocation des Bords Tampons  **************************/
    /***************************************************************************/
    DimensionMax=max(X,Y) + 1;
    Carte_Tampon= (PREC_CARTE *)  calloc(DimensionMax ,sizeof(PREC_CARTE));
    BufferBord01= (PREC_GRILLE *) calloc(DimensionMax,sizeof(PREC_GRILLE));
    BufferBord02= (PREC_GRILLE *) calloc(W+1,sizeof(PREC_GRILLE));
    BufferBord03= (PREC_GRILLE *) calloc(W+1,sizeof(PREC_GRILLE));
    /***************************************************************************/
    /**************** Allocation des Bords PML *********************************/
    /***************************************************************************/
    if(ConditionBordXp==0)
    {
        Vx_Xp_X =(PREC_GRILLE *) calloc(W*Y,sizeof(PREC_GRILLE));
        Vx_Xp_Y =(PREC_GRILLE *) calloc(W*Y,sizeof(PREC_GRILLE));
        Vy_Xp_X =(PREC_GRILLE *) calloc((W+1)*(Y+1),sizeof(PREC_GRILLE));
        Vy_Xp_Y =(PREC_GRILLE *) calloc((W+1)*(Y+1),sizeof(PREC_GRILLE));
        Txx_Xp_X=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Txx_Xp_Y=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Tyy_Xp_X=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Tyy_Xp_Y=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Txy_Xp_X=(PREC_GRILLE *) calloc(W*(Y+1),sizeof(PREC_GRILLE));
        Txy_Xp_Y=(PREC_GRILLE *) calloc(W*(Y+1),sizeof(PREC_GRILLE));
    }
    if(ConditionBordXm==0)
    {
        Vx_Xm_X =(PREC_GRILLE *) calloc(W*Y,sizeof(PREC_GRILLE));
        Vx_Xm_Y =(PREC_GRILLE *) calloc(W*Y,sizeof(PREC_GRILLE));
        Vy_Xm_X =(PREC_GRILLE *) calloc((W+1)*(Y+1),sizeof(PREC_GRILLE));
        Vy_Xm_Y =(PREC_GRILLE *) calloc((W+1)*(Y+1),sizeof(PREC_GRILLE));
        Txx_Xm_X=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Txx_Xm_Y=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Tyy_Xm_X=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Tyy_Xm_Y=(PREC_GRILLE *) calloc((W+1)*Y,sizeof(PREC_GRILLE));
        Txy_Xm_X=(PREC_GRILLE *) calloc(W*(Y+1),sizeof(PREC_GRILLE));
        Txy_Xm_Y=(PREC_GRILLE *) calloc(W*(Y+1),sizeof(PREC_GRILLE));
    }
    if(ConditionBordYp==0)
    {
        Vy_Yp_X =(PREC_GRILLE *) calloc(W*X,sizeof(PREC_GRILLE));
        Vy_Yp_Y =(PREC_GRILLE *) calloc(W*X,sizeof(PREC_GRILLE));
        Vx_Yp_X =(PREC_GRILLE *) calloc((W+1)*(X+1),sizeof(PREC_GRILLE));
        Vx_Yp_Y =(PREC_GRILLE *) calloc((W+1)*(X+1),sizeof(PREC_GRILLE));
        Tyy_Yp_X=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Tyy_Yp_Y=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Txx_Yp_X=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Txx_Yp_Y=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Txy_Yp_X=(PREC_GRILLE *) calloc(W*(X+1),sizeof(PREC_GRILLE));
        Txy_Yp_Y=(PREC_GRILLE *) calloc(W*(X+1),sizeof(PREC_GRILLE));
    }
    if(ConditionBordYm==0)
    {
        Vy_Ym_X =(PREC_GRILLE *) calloc(W*X,sizeof(PREC_GRILLE));
        Vy_Ym_Y =(PREC_GRILLE *) calloc(W*X,sizeof(PREC_GRILLE));
        Vx_Ym_X =(PREC_GRILLE *) calloc((W+1)*(X+1),sizeof(PREC_GRILLE));
        Vx_Ym_Y =(PREC_GRILLE *) calloc((W+1)*(X+1),sizeof(PREC_GRILLE));
        Txx_Ym_X=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Txx_Ym_Y=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Tyy_Ym_X=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Tyy_Ym_Y=(PREC_GRILLE *) calloc((W+1)*X,sizeof(PREC_GRILLE));
        Txy_Ym_X=(PREC_GRILLE *) calloc(W*(X+1),sizeof(PREC_GRILLE));
        Txy_Ym_Y=(PREC_GRILLE *) calloc(W*(X+1),sizeof(PREC_GRILLE));
    }
    /***************************************************************************/
    /**************** Allocation des Coins PML *********************************/
    /***************************************************************************/
    if((ConditionBordXp==0) && (ConditionBordYp==0))
    {
        Vx_XpYp_X =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vx_XpYp_Y =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vy_XpYp_X =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Vy_XpYp_Y =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Txx_XpYp_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txx_XpYp_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XpYp_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XpYp_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txy_XpYp_X=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
        Txy_XpYp_Y=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
    }
    if((ConditionBordXp==0) && (ConditionBordYm==0))
    {
        Vx_XpYm_X =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vx_XpYm_Y =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vy_XpYm_X =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Vy_XpYm_Y =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Txx_XpYm_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txx_XpYm_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XpYm_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XpYm_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txy_XpYm_X=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
        Txy_XpYm_Y=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
    }
    if((ConditionBordXm==0) && (ConditionBordYp==0))
    {
        Vx_XmYp_X =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vx_XmYp_Y =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vy_XmYp_X =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Vy_XmYp_Y =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Txx_XmYp_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txx_XmYp_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XmYp_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XmYp_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txy_XmYp_X=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
        Txy_XmYp_Y=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
    }
    if((ConditionBordXm==0) && (ConditionBordYm==0))
    {
        Vx_XmYm_X =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vx_XmYm_Y =(PREC_GRILLE *) calloc((W)*(W+1),sizeof(PREC_GRILLE));
        Vy_XmYm_X =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Vy_XmYm_Y =(PREC_GRILLE *) calloc((W+1)*(W),sizeof(PREC_GRILLE));
        Txx_XmYm_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txx_XmYm_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XmYm_X=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Tyy_XmYm_Y=(PREC_GRILLE *) calloc((W+1)*(W+1),sizeof(PREC_GRILLE));
        Txy_XmYm_X=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
        Txy_XmYm_Y=(PREC_GRILLE *) calloc((W)*(W),sizeof(PREC_GRILLE));
    }
}
void LiberationMemoire(void)
{

    free(Carte);
    free(AttenuationPML);

    /***************************************************************************/
    /*********** Liberation des grilles du bloc principal **********************/
    /***************************************************************************/
    free(Vx);
    free(Vy);
    free(Txx);
    free(Tyy);
    free(Txy);
    if(TypeAbsorption==1)
    {
        free(Vx_Old);
        free(Vy_Old);
    }
    if(SnapshotUx==1)
        free(Ux);
    if(SnapshotUy==1)
        free(Uy);
    if(SnapshotMaxTxx==1)
        free(MaxTxx);
    if(SnapshotMaxTyy==1)
        free(MaxTyy);
    if(SnapshotMaxTxy==1)
        free(MaxTxy);
    if(SnapshotMaxVx==1)
        free(MaxVx);
    if(SnapshotMaxVy==1)
        free(MaxVy);
    if(SnapshotMaxV==1)
        free(MaxV);
    if(SnapshotMaxEc==1)
        free(MaxEc);
    if(SnapshotMaxCurl==1)
        free(MaxCurl);
    if(SnapshotMaxDiv==1)
        free(MaxDiv);
    /***************************************************************************/
    /**************** Liberation des Bords PML *********************************/
    /***************************************************************************/
    if(ConditionBordXp==0)
    {
        free(Vx_Xp_X);
        free(Vx_Xp_Y);
        free(Vy_Xp_X);
        free(Vy_Xp_Y);
        free(Txx_Xp_X);
        free(Txx_Xp_Y);
        free(Tyy_Xp_X);
        free(Tyy_Xp_Y);
        free(Txy_Xp_X);
        free(Txy_Xp_Y);
    }
    if(ConditionBordXm==0)
    {
        free(Vx_Xm_X);
        free(Vx_Xm_Y);
        free(Vy_Xm_X);
        free(Vy_Xm_Y);
        free(Txx_Xm_X);
        free(Txx_Xm_Y);
        free(Tyy_Xm_X);
        free(Tyy_Xm_Y);
        free(Txy_Xm_X);
        free(Txy_Xm_Y);
    }
    if(ConditionBordYm==0)
    {
        free(Vy_Ym_X);
        free(Vy_Ym_Y);
        free(Vx_Ym_X);
        free(Vx_Ym_Y);
        free(Txx_Ym_X);
        free(Txx_Ym_Y);
        free(Tyy_Ym_X);
        free(Tyy_Ym_Y);
        free(Txy_Ym_X);
        free(Txy_Ym_Y);
    }
    if(ConditionBordYp==0)
    {
        free(Vy_Yp_X);
        free(Vy_Yp_Y);
        free(Txy_Yp_X);
        free(Txy_Yp_Y);
        free(Txx_Yp_X);
        free(Txx_Yp_Y);
        free(Vx_Yp_X);
        free(Vx_Yp_Y);
        free(Tyy_Yp_X);
        free(Tyy_Yp_Y);
    }

    /***************************************************************************/
    /**************** Liberation des Coins PML *********************************/
    /***************************************************************************/
    if((ConditionBordXp==0) && (ConditionBordYp==0))
    {
        free(Vx_XpYp_X);
        free(Vx_XpYp_Y);
        free(Vy_XpYp_X);
        free(Vy_XpYp_Y);
        free(Txx_XpYp_X);
        free(Txx_XpYp_Y);
        free(Tyy_XpYp_X);
        free(Tyy_XpYp_Y);
        free(Txy_XpYp_X);
        free(Txy_XpYp_Y);
    }
    if((ConditionBordXp==0) && (ConditionBordYm==0))
    {
        free(Vx_XpYm_X);
        free(Vx_XpYm_Y);
        free(Vy_XpYm_X);
        free(Vy_XpYm_Y);
        free(Txx_XpYm_X);
        free(Txx_XpYm_Y);
        free(Tyy_XpYm_X);
        free(Tyy_XpYm_Y);
        free(Txy_XpYm_X);
        free(Txy_XpYm_Y);
    }
    if((ConditionBordXm==0) && (ConditionBordYp==0))
    {
        free(Vx_XmYp_X);
        free(Vx_XmYp_Y);
        free(Vy_XmYp_X);
        free(Vy_XmYp_Y);
        free(Txx_XmYp_X);
        free(Txx_XmYp_Y);
        free(Tyy_XmYp_X);
        free(Tyy_XmYp_Y);
        free(Txy_XmYp_X);
        free(Txy_XmYp_Y);
    }
    if((ConditionBordXm==0) && (ConditionBordYm==0))
    {
        free(Vx_XmYm_X);
        free(Vx_XmYm_Y);
        free(Vy_XmYm_X);
        free(Vy_XmYm_Y);
        free(Txx_XmYm_X);
        free(Txx_XmYm_Y);
        free(Tyy_XmYm_X);
        free(Tyy_XmYm_Y);
        free(Txy_XmYm_X);
        free(Txy_XmYm_Y);
    }
    /***************************************************************************/
    /****************** Liberation des Bords Tampons  **************************/
    /***************************************************************************/
    free(Carte_Tampon);
    free(BufferBord01);
    free(BufferBord02);
    free(BufferBord03);
    free(TableLegerete);
    free(TableC11);
    free(TableC22);
    free(TableC12);
    free(TableC33);
    free(TableX11);
    free(TableX22);
    free(TableX12);
    free(TableX33);
    free(EmissionTxx);
    free(EmissionTyy);
    free(EmissionTxy);
    free(EmissionVx);
    free(EmissionVy);
    /*
    A faire proprement la libération de mémoire...

     free(SignauxReceptionTxx);
     free(SignauxReceptionTyy);
     free(SignauxReceptionTxy);
     free(SignauxReceptionVx);
     free(SignauxReceptionVy);
    */
    free(ReceptionTxx);
    free(ReceptionTyy);
    free(ReceptionTxy);
    free(ReceptionVx);
    free(ReceptionVy);
}

/************************ Affichage  ***********************************************/

void AffichageTempsRestant(void)
{
    int Heures, Minutes, Secondes;
#ifdef WIN32
    int i;
#endif

    TempsCourant=time(NULL);
    TempsRestant=difftime(TempsCourant,TempsDebut)*((double)(IndiceTempsMax+1)/(double)(IndiceTemps+1)-1);
    Heures=(int)(TempsRestant/3600);
    Minutes=(int)((TempsRestant-3600*Heures)/60);
    Secondes=(int)(TempsRestant-3600*Heures-60*Minutes);
#ifdef WIN32
    for(i=0; i < 80; i++)
        printf("\b");
#else
    printf("\n");
#endif

    printf("Computed: %.1f/%.1f us <--> ",IndiceTemps*DeltaT,(IndiceTempsMax-1)*DeltaT);
    printf("Step: %i/%i   ",IndiceTemps,IndiceTempsMax-1);
    printf("Remaining time: %ih %im %is ",Heures, Minutes,Secondes);
}
void AffichageDureeSimul(void)
{
    int Heures, Minutes, Secondes;
    TempsCourant=time(NULL);
    TempsRestant=difftime(TempsCourant,TempsDebut);
    Heures=(int)(TempsRestant/3600);
    Minutes=(int)((TempsRestant-3600*Heures)/60);
    Secondes=(int)(TempsRestant-3600*Heures-60*Minutes);
    printf("Total computation time: %3ih %2imin %2isec \n",Heures, Minutes,Secondes);
}

void AffichageStartTime(void)
{

    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf ( "Started on : %s", asctime (timeinfo) );
}

void AffichageEndTime(void)
{

    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf ( "\nEnded on   : %s", asctime (timeinfo) );
}


void MAJBarretteRecepteurs(void)
{
    int NumeroBarrette;
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionTxx; NumeroBarrette++)
    {
        EnregistrementBarretteReception(Txx,&ReceptionTxx[NumeroBarrette],X,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionTyy; NumeroBarrette++)
    {
        EnregistrementBarretteReception(Tyy,&ReceptionTyy[NumeroBarrette],X,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionTxy; NumeroBarrette++)
    {
        EnregistrementBarretteReception(Txy,&ReceptionTxy[NumeroBarrette],X+1,Y+1);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionVx; NumeroBarrette++)
    {
        EnregistrementBarretteReception(Vx,&ReceptionVx[NumeroBarrette],X+1,Y);
    }
    for(NumeroBarrette=0; NumeroBarrette < NbBarrettesReceptionVy; NumeroBarrette++)
    {
        EnregistrementBarretteReception(Vy,&ReceptionVy[NumeroBarrette],X,Y+1);
    }
}

/*********************************** Manipulation de SIGNAUX **********************/

void GenerationSignal(int TypeSignal, double Freq, double T0, double Bandwidth,double *Signal)
{
    int i,i0;
    double argument,Maximum;
    Maximum=0.0;
    i0=(int)(T0/DeltaT);
    switch(TypeSignal)
    {
    case 1:
    default:
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=2*Bandwidth*Freq*(double)(i-i0)*DeltaT;
            Signal[i]=sin(2*3.1415*Freq*(double)(i-i0)*DeltaT)*exp(-argument*argument);
        }
        break;
    case 2 : 	/*Signal à moyenne nulle, de primitive à moyenne non nulle*/
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=Freq*(double)(i-i0)*DeltaT;
            Signal[i]=- (i-i0)*exp(-argument*argument);
        }
        break;
    case 3 :
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=2*Bandwidth*Freq*(double)(i-i0)*DeltaT;
            Signal[i]=cos(2*3.1415*Freq*(double)(i-i0)*DeltaT)*exp(-argument*argument);
        }
        break;
    case 10 : 	/*Signal à moyenne nulle, de primitive à moyenne non nulle*/
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=4.44*Freq*(double)(i-i0)*DeltaT;
            Signal[i]=- (i-i0)*exp(-argument*argument);
        }
        break;
    case 11 :	/*Dérivé du précedent*/
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=3.116*Freq*(double)(i-i0)*DeltaT;
            Signal[i]=- (i-i0)*exp(-argument*argument);
        }
        for(i=0; i < (IndiceTempsMax-1); i++)  /* dérivation*/
        {
            Signal[i]=Signal[i+1]-Signal[i];
        }
        break;
    case 12 :	/*Dérivé du précedent*/
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=2.5*Freq*(double)(i-i0)*DeltaT;
            Signal[i]=- (i-i0)*exp(-argument*argument);
        }
        for(i=0; i < (IndiceTempsMax-2); i++)  /* dérivation deux fois...*/
        {
            Signal[i]=Signal[i+2]+Signal[i]-2*Signal[i+1];
        }
        break;
    case 13 : 	/*Dérivé du précedent*/
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=2.2*Freq*(double)(i-i0)*DeltaT;
            Signal[i]=- (i-i0)*exp(-argument*argument);
        }
        for(i=0; i < (IndiceTempsMax-3); i++)   /*dérivation trois fois...*/
        {
            Signal[i]=Signal[i+3]-3*Signal[i+2]+3*Signal[i+1]-Signal[i];
        }
        break;
    case 14 :	/*Dérivé du précedent*/
        for(i=0; i < IndiceTempsMax; i++)
        {
            argument=2.05*Freq*(double)(i-i0)*DeltaT;
            Signal[i]=-(i-i0)*exp(-argument*argument);
        }
        for(i=0; i < (IndiceTempsMax-4); i++)  /* dérivation trois fois...*/
        {
            Signal[i]=Signal[i+4]-4*Signal[i+3]+6*Signal[i+2]-4*Signal[i+1]+Signal[i];
        }
        break;
    }
    /*Normalisation du signal*/
    Maximum=MaxAbsSignal(Signal,IndiceTempsMax);
    if(!(Maximum==0))
        for(i=0; i < IndiceTempsMax; i++)
        {

            Signal[i]=Signal[i]/Maximum;
        }
    /**/
}
double MaxAbsSignal(double* Signal,int NbptsSignal)
{
    int i;
    double Maximum=0;
    for(i=0; i<NbptsSignal; i++)
        Maximum=max(Maximum,fabs(Signal[i]));
    return Maximum;
}

/******************************GESTION des TAMPONS *********************************/

void MiseAJourContraintesTamponBordPML(const char *TypeBord)
{
    int i,j;
    if(!strcmp(TypeBord,"Xp"))
    {
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y; j++)
        {
            Carte_Tampon[j]=Carte[(X-1)*Y+j];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y+1; j++)
        {
            Tij_BufferBordBlocI[j]=Txy[X*(Y+1)+j];
        }
        if((ConditionBordYp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJp[i]=Tyy_XpYp_X[(i)*(W+1)]	+ Tyy_XpYp_Y[(i)*(W+1)];
            }
        }
        if((ConditionBordYm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJm[i]=Tyy_XpYm_X[(i)*(W+1)]	+ Tyy_XpYm_Y[(i)*(W+1)];
            }
        }
    }
    if(!strcmp(TypeBord,"Xm"))
    {
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y; j++)
        {
            Carte_Tampon[j]=Carte[0*Y+j];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y+1; j++)
        {
            Tij_BufferBordBlocI[j]=Txy[0*(Y+1)+j];
        }
        if((ConditionBordYp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJp[i]=Tyy_XmYp_X[(i)*(W+1)]	+ Tyy_XmYp_Y[(i)*(W+1)];
            }
        }
        if((ConditionBordYm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJm[i]=Tyy_XmYm_X[(i)*(W+1)]	+ Tyy_XmYm_Y[(i)*(W+1)];
            }
        }
    }
    if(!strcmp(TypeBord,"Yp"))
    {
#pragma omp parallel for schedule(guided)
        for(j=0; j < X; j++)
        {
            Carte_Tampon[j]=Carte[j*Y+Y-1];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < X+1; j++)
        {
            Tij_BufferBordBlocI[j]=Txy[j*(Y+1)+Y];
        }
        if((ConditionBordXp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJp[i]=Txx_XpYp_X[i]	+ Txx_XpYp_Y[i];
            }
        }
        if((ConditionBordXm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJm[i]=Txx_XmYp_X[i]	+ Txx_XmYp_Y[i];
            }
        }
    }
    if(!strcmp(TypeBord,"Ym"))
    {
#pragma omp parallel for schedule(guided)
        for(j=0; j < X; j++)
        {
            Carte_Tampon[j]=Carte[j*Y+0];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < X+1; j++)
        {
            Tij_BufferBordBlocI[j]=Txy[j*(Y+1)+0];
        }
        if((ConditionBordXp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJp[i]=Txx_XpYm_X[i]	+ Txx_XpYm_Y[i];
            }
        }
        if((ConditionBordXm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W+1; i++)
            {
                Tjj_BufferBordCoinJm[i]=Txx_XmYm_X[i]	+ Txx_XmYm_Y[i];
            }
        }
    }
}
void MiseAJourVitessesTamponBordPML(const char *TypeBord)
{
    int i,j;
    if(!strcmp(TypeBord,"Xp"))
    {
        TableCii=TableC11;
        TableCjj=TableC22;
        TableCij=TableC12;
        TableCijij=TableC33;
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y; j++)
        {
            Carte_Tampon[j]=Carte[(X-1)*Y+j];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y; j++)
        {
            Vi_BufferBordBlocI[j]=Vx[(X)*(Y)+j];
        }
        if((ConditionBordYp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJp[i]=Vx_XpYp_X[(i)*(W+1)] +	Vx_XpYp_Y[(i)*(W+1)];
            }
        }
        if((ConditionBordYm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJm[i]=Vx_XpYm_X[(i)*(W+1)] +	Vx_XpYm_Y[(i)*(W+1)];
            }
        }
    }
    if(!strcmp(TypeBord,"Xm"))
    {
        TableCii=TableC11;
        TableCjj=TableC22;
        TableCij=TableC12;
        TableCijij=TableC33;
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y; j++)
        {
            Carte_Tampon[j]=Carte[(0)*Y+j];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < Y; j++)
        {
            Vi_BufferBordBlocI[j]=Vx[(0)*(Y)+j];
        }
        if((ConditionBordYp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJp[i]=Vx_XmYp_X[(i)*(W+1)] +	Vx_XmYp_Y[(i)*(W+1)];
            }
        }
        if((ConditionBordYm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJm[i]=Vx_XmYm_X[(i)*(W+1)] +	Vx_XmYm_Y[(i)*(W+1)];
            }
        }
    }
    if(!strcmp(TypeBord,"Yp"))
    {
        TableCii=TableC22;
        TableCjj=TableC11;
        TableCij=TableC12;
        TableCijij=TableC33;
#pragma omp parallel for schedule(guided)
        for(j=0; j < X; j++)
        {
            Carte_Tampon[j]=Carte[(j)*Y+Y-1];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < X; j++)
        {
            Vi_BufferBordBlocI[j]=Vy[(j)*(Y+1)+Y];
        }
        if((ConditionBordXp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJp[i]=Vy_XpYp_X[i] +	Vy_XpYp_Y[i];
            }
        }
        if((ConditionBordXm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJm[i]=Vy_XmYp_X[i] +	Vy_XmYp_Y[i];
            }
        }
    }
    if(!strcmp(TypeBord,"Ym"))
    {
        TableCii=TableC22;
        TableCjj=TableC11;
        TableCij=TableC12;
        TableCijij=TableC33;
#pragma omp parallel for schedule(guided)
        for(j=0; j < X; j++)
        {
            Carte_Tampon[j]=Carte[(j)*Y+0];
        }
#pragma omp parallel for schedule(guided)
        for(j=0; j < X; j++)
        {
            Vi_BufferBordBlocI[j]=Vy[(j)*(Y+1)+0];
        }
        if((ConditionBordXp==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJp[i]=Vy_XpYm_X[i] +	Vy_XpYm_Y[i];
            }
        }
        if((ConditionBordXm==0))
        {
#pragma omp parallel for schedule(guided)
            for(i=0; i < W; i++)
            {
                Vi_BufferBordCoinJm[i]=Vy_XmYm_X[i] +	Vy_XmYm_Y[i];
            }
        }
    }
}
void MiseAJourContraintesTamponCoinPML(const char *TypeCoin)
{
    int i;
    if(!strcmp(TypeCoin,"XpYp"))
    {
        Carte_Tampon[0]=Carte[(X-1)*Y+Y-1];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordI[i]=Txy_Xp_X[(i)*(Y+1)+Y] + Txy_Xp_Y[(i)*(Y+1)+Y];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordJ[i]=Txy_Yp_X[(i)*(X+1)+X] + Txy_Yp_Y[(i)*(X+1)+X];
        }
    }
    if(!strcmp(TypeCoin,"XmYp"))
    {
        Carte_Tampon[0]=Carte[(0)*Y+Y-1];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordI[i]=Txy_Xm_X[(i)*(Y+1)+Y] + Txy_Xm_Y[(i)*(Y+1)+Y];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordJ[i]=Txy_Yp_X[(i)*(X+1)+0] + Txy_Yp_Y[(i)*(X+1)+0];
        }
    }
    if(!strcmp(TypeCoin,"XpYm"))
    {
        Carte_Tampon[0]=Carte[(X-1)*Y+0];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordI[i]=Txy_Xp_X[(i)*(Y+1)+0] + Txy_Xp_Y[(i)*(Y+1)+0];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordJ[i]=Txy_Ym_X[(i)*(X+1)+X] + Txy_Ym_Y[(i)*(X+1)+X];
        }
    }
    if(!strcmp(TypeCoin,"XmYm"))
    {
        Carte_Tampon[0]=Carte[(0)*Y+0];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordI[i]=Txy_Xm_X[(i)*(Y+1)+0] + Txy_Xm_Y[(i)*(Y+1)+0];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W; i++)
        {
            Tij_BufferCoinBordJ[i]=Txy_Ym_X[(i)*(X+1)+0] + Txy_Ym_Y[(i)*(X+1)+0];
        }
    }
}
void MiseAJourVitessesTamponCoinPML(const char *TypeCoin)
{
    int i;
    /*pas de rotation dans les coins...*/
    TableCii=TableC11;
    TableCjj=TableC22;
    TableCij=TableC12;
    TableCijij=TableC33;
    if(!strcmp(TypeCoin,"XpYp"))
    {
        Carte_Tampon[0]=Carte[(X-1)*Y+Y-1];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vj_BufferCoinBordI[i]=	Vy_Xp_X[(i)*(Y+1)+Y] + Vy_Xp_Y[(i)*(Y+1)+Y];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vi_BufferCoinBordJ[i]=	Vx_Yp_X[(i)*(X+1)+X] + Vx_Yp_Y[(i)*(X+1)+X];
        }
    }
    if(!strcmp(TypeCoin,"XpYm"))
    {
        Carte_Tampon[0]=Carte[(X-1)*Y+0];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vj_BufferCoinBordI[i]=	Vy_Xp_X[(i)*(Y+1)+0] + Vy_Xp_Y[(i)*(Y+1)+0];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vi_BufferCoinBordJ[i]=	Vx_Ym_X[(i)*(X+1)+X] + Vx_Ym_Y[(i)*(X+1)+X];
        }
    }
    if(!strcmp(TypeCoin,"XmYp"))
    {
        Carte_Tampon[0]=Carte[(0)*Y+Y-1];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vj_BufferCoinBordI[i]=	Vy_Xm_X[(i)*(Y+1)+Y] + Vy_Xm_Y[(i)*(Y+1)+Y];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vi_BufferCoinBordJ[i]=	Vx_Yp_X[(i)*(X+1)+0] + Vx_Yp_Y[(i)*(X+1)+0];
        }
    }
    if(!strcmp(TypeCoin,"XmYm"))
    {
        Carte_Tampon[0]=Carte[(0)*Y+0];
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vj_BufferCoinBordI[i]=	Vy_Xm_X[(i)*(Y+1)+0] + Vy_Xm_Y[(i)*(Y+1)+0];
        }
#pragma omp parallel for schedule(guided)
        for(i=0; i < W+1; i++)
        {
            Vi_BufferCoinBordJ[i]=	Vx_Ym_X[(i)*(X+1)+0] + Vx_Ym_Y[(i)*(X+1)+0];
        }
    }
}
void ReMiseAZeroTampons(void)
{
    int DimensionMax,i;
    DimensionMax=max(X,Y)+1;
#pragma omp parallel for schedule(guided)
    for(i=0; i < DimensionMax; i++)
    {
        Carte_Tampon[i]=0;
        BufferBord01[i]=0;
    }
#pragma omp parallel for schedule(guided)
    for(i=0; i < W+1; i++)
    {
        BufferBord02[i]=0;
        BufferBord03[i]=0;
    }
}
void AdressageBuffer(void)
{
    Vi_BufferBordBlocI		=BufferBord01;
    Vi_BufferBordCoinJp 	=BufferBord02;
    Vi_BufferBordCoinJm 	=BufferBord03;
    Tij_BufferBordBlocI 	=BufferBord01;
    Tjj_BufferBordCoinJp 	=BufferBord02;
    Tjj_BufferBordCoinJm 	=BufferBord03;
    Vj_BufferCoinBordI		=BufferBord02;
    Vi_BufferCoinBordJ		=BufferBord03;
    Tij_BufferCoinBordI 	=BufferBord02;
    Tij_BufferCoinBordJ 	=BufferBord03;
}
