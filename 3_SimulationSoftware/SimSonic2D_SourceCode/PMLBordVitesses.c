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

/** \file PMLBordVitesses.c
    \brief Code source des fonctions de calcul des PML pour les vitesses sur les bords de la carte.
    \author Emmanuel Bossy
    \version 2D_2011_10_28_omp
    \date 28 octobre 2011 */

void CalculPMLBordVitesses(PREC_GRILLE *Vi_I,PREC_GRILLE *Vi_J,
                           PREC_GRILLE *Vj_I,PREC_GRILLE *Vj_J, PREC_GRILLE *Tii_I,PREC_GRILLE *Tii_J,
                           PREC_GRILLE *Tjj_I,PREC_GRILLE *Tjj_J, PREC_GRILLE *Tij_I,PREC_GRILLE *Tij_J,
                           int J, const char *TypeBord)
{
    int i,j ;
    double Legerete ;
    double FacteurQModel ;
    double SignI ;
    if(!strcmp(TypeBord,"Xp") || !strcmp(TypeBord,"Yp"))
        SignI = +1.0 ;
    else
        SignI = -1.0 ;
    /*************************************************************************/
    /*****************************Bloc du Bord PML****************************/
    /*************************************************************************/
	#pragma omp parallel for private(j,Legerete) schedule(guided)
    for(i = 0 ; i < W ; i++)
    {
        for(j = 0 ; j < J ; j++)
        {
            Legerete = TableLegerete[Carte_Tampon[j]] ;
            Vi_I[i*(J)+j] = 1./(1+DeltaT*AttenuationPML[1+2*i])*(
                                (1-DeltaT*AttenuationPML[1+2*i])*Vi_I[i*(J)+j]+ DeltaTsurH*Legerete*SignI*(
                                    (Tii_I[(i+1)*(J)+j] + Tii_J[(i+1)*(J)+j]) - (Tii_I[(i)*(J)+j] + Tii_J[(i)*(J)+j]))) ;
            Vi_J[i*(J)+j] += DeltaTsurH*Legerete*(
                                 (Tij_I[i*(J+1)+j+1] + Tij_J[i*(J+1)+j+1]) - (Tij_I[i*(J+1)+j] + Tij_J[i*(J+1)+j])) ;
        }
    }
	#pragma omp parallel for private(j,Legerete) schedule(guided)
    for(i = 1 ; i < W ; i++)
    {
        for(j = 1 ; j < J ; j++)
        {
            Legerete = 2./(1./TableLegerete[Carte_Tampon[j-1]]+1./TableLegerete[Carte_Tampon[j]]) ;
            Vj_I[i*(J+1)+j] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                                  (1-DeltaT*AttenuationPML[2*i])*Vj_I[i*(J+1)+j]+ DeltaTsurH*Legerete*SignI*(
                                      (Tij_I[(i)*(J+1)+j] + Tij_J[(i)*(J+1)+j]) - (Tij_I[(i-1)*(J+1)+j] + Tij_J[(i-1)*(J+1)+j]))) ;
            Vj_J[i*(J+1)+j] += DeltaTsurH*Legerete*(
                                   (Tjj_I[i*J+j] + Tjj_J[i*J+j]) - (Tjj_I[i*J+j-1] + Tjj_J[i*J+j-1])) ;
        }
    }
    /*************************************************************************/
    /****************************Bord sur Bloc du Bord PML ...****************/
    /*************************************************************************/
	#pragma omp parallel for private(Legerete) schedule(guided)
    for(j = 1 ; j < J ; j++)
    {

        Legerete = 2./(1./TableLegerete[Carte_Tampon[j-1]]+1./TableLegerete[Carte_Tampon[j]]) ;

        Vj_I[0*(J+1)+j] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                              (1-DeltaT*AttenuationPML[2*0])*Vj_I[0*(J+1)+j]+ DeltaTsurH*Legerete*SignI*(
                                  (Tij_I[(0)*(J+1)+j] + Tij_J[(0)*(J+1)+j]) - (Tij_BufferBordBlocI[j]))) ;
        Vj_J[0*(J+1)+j] += DeltaTsurH*Legerete*(
                               (Tjj_I[0*J+j] +	Tjj_J[0*J+j]) - (Tjj_I[0*J+j-1]	+ Tjj_J[0*J+j-1])) ;
    }
    /*************************************************************************/
    /****************************Bord sur Coin du Bord PML ...****************/
    /*************************************************************************/
    /****************************Coin PML ...*********************************/
    if(((!strcmp(TypeBord,"Xp")) && (ConditionBordYp == 0)) ||
            ((!strcmp(TypeBord,"Xm")) && (ConditionBordYp == 0)) ||
            ((!strcmp(TypeBord,"Yp")) && (ConditionBordXp == 0)) ||
            ((!strcmp(TypeBord,"Ym")) && (ConditionBordXp == 0)))
    {
		#pragma omp parallel for private(Legerete) schedule(guided)
        for(i = 1 ; i < W ; i++)
        {
            Legerete = (TableLegerete[Carte_Tampon[J-1]]) ;
            Vj_I[i*(J+1)+J] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                                  (1-DeltaT*AttenuationPML[2*i])*Vj_I[i*(J+1)+J]+ DeltaTsurH*Legerete*SignI*(
                                      (Tij_I[(i)*(J+1)+J] + Tij_J[(i)*(J+1)+J]) - (Tij_I[(i-1)*(J+1)+J] + Tij_J[(i-1)*(J+1)+J]))) ;
            Vj_J[i*(J+1)+J] += DeltaTsurH*Legerete*(
                                   (Tjj_BufferBordCoinJp[i]) - (Tjj_I[i*J+J-1] + Tjj_J[i*J+J-1])) ;
        }
        /*Coin du bord PML*/
        Legerete = (TableLegerete[Carte_Tampon[J-1]]) ;
        Vj_I[0*(J+1)+J] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                              (1-DeltaT*AttenuationPML[2*0])*Vj_I[0*(J+1)+J]+ DeltaTsurH*Legerete*SignI*(
                                  (Tij_I[(0)*(J+1)+J] + Tij_J[(0)*(J+1)+J]) - (Tij_BufferBordBlocI[J]))) ;
        Vj_J[0*(J+1)+J] += DeltaTsurH*Legerete*(
                               (Tjj_BufferBordCoinJp[0]) - (Tjj_I[0*J+J-1] + Tjj_J[0*J+J-1])) ;
    }
    if(((!strcmp(TypeBord,"Xp")) && (ConditionBordYm == 0)) ||
            ((!strcmp(TypeBord,"Xm")) && (ConditionBordYm == 0)) ||
            ((!strcmp(TypeBord,"Yp")) && (ConditionBordXm == 0)) ||
            ((!strcmp(TypeBord,"Ym")) && (ConditionBordXm == 0)))
    {
		#pragma omp parallel for private(Legerete) schedule(guided)
        for(i = 1 ; i < W ; i++)
        {
            Legerete = (TableLegerete[Carte_Tampon[0]]) ;
            Vj_I[i*(J+1)+0] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                                  (1-DeltaT*AttenuationPML[2*i])*Vj_I[i*(J+1)+0]+ DeltaTsurH*Legerete*SignI*(
                                      (Tij_I[(i)*(J+1)+0] + Tij_J[(i)*(J+1)+0]) - (Tij_I[(i-1)*(J+1)+0] + Tij_J[(i-1)*(J+1)+0]))) ;
            Vj_J[i*(J+1)+0] += DeltaTsurH*Legerete*(
                                   (Tjj_I[i*J+0] +	Tjj_J[i*J+0]) -	(Tjj_BufferBordCoinJm[i])) ;
        }
        /*Coin du bord PML*/
        Legerete = (TableLegerete[Carte_Tampon[0]]) ;
        Vj_I[0*(J+1)+0] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                              (1-DeltaT*AttenuationPML[2*0])*Vj_I[0*(J+1)+0]+ DeltaTsurH*Legerete*SignI*(
                                  (Tij_I[(0)*(J+1)+0] + Tij_J[(0)*(J+1)+0]) - (Tij_BufferBordBlocI[0]))) ;
        Vj_J[0*(J+1)+0] += DeltaTsurH*Legerete*(
                               (Tjj_I[0*J+0] +	Tjj_J[0*J+0]) -	(Tjj_BufferBordCoinJm[0])) ;
    }

    /*****************Frontières Coin Libre...*********************************/
    if(((!strcmp(TypeBord,"Xp")) && ((ConditionBordYp == 2)||(ConditionBordYp == 4))) ||
            ((!strcmp(TypeBord,"Xm")) && ((ConditionBordYp == 2)||(ConditionBordYp == 4))) ||
            ((!strcmp(TypeBord,"Yp")) && ((ConditionBordXp == 2)||(ConditionBordXp == 4))) ||
            ((!strcmp(TypeBord,"Ym")) && ((ConditionBordXp == 2)||(ConditionBordXp == 4))))
    {
		#pragma omp parallel for private(Legerete) schedule(guided)
        for(i = 1 ; i < W ; i++)
        {
            Legerete = (TableLegerete[Carte_Tampon[J-1]]);
            Vj_I[i*(J+1)+J] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                                  (1-DeltaT*AttenuationPML[2*i])*Vj_I[i*(J+1)+J]+ DeltaTsurH*Legerete*SignI*(
                                      (Tij_I[(i)*(J+1)+J] + Tij_J[(i)*(J+1)+J]) - (Tij_I[(i-1)*(J+1)+J] + Tij_J[(i-1)*(J+1)+J]))) ;
            Vj_J[i*(J+1)+J] += DeltaTsurH*Legerete*(
                                   0 - 2*(Tjj_I[i*J+J-1] +	Tjj_J[i*J+J-1])) ;
        }
        /*Coin du bord PML*/
        Legerete = (TableLegerete[Carte_Tampon[J-1]]);
        Vj_I[0*(J+1)+J] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                              (1-DeltaT*AttenuationPML[2*0])*Vj_I[0*(J+1)+J]+ DeltaTsurH*Legerete*SignI*(
                                  (Tij_I[(0)*(J+1)+J] + Tij_J[(0)*(J+1)+J]) - (Tij_BufferBordBlocI[J]))) ;
        Vj_J[0*(J+1)+J] += DeltaTsurH*Legerete*(
                               0 - 2*(Tjj_I[0*J+J-1] +	Tjj_J[0*J+J-1])) ;
    }
    if(((!strcmp(TypeBord,"Xp")) && ((ConditionBordYm == 2)||(ConditionBordYm == 4))) ||
            ((!strcmp(TypeBord,"Xm")) && ((ConditionBordYm == 2)||(ConditionBordYm == 4))) ||
            ((!strcmp(TypeBord,"Yp")) && ((ConditionBordXm == 2)||(ConditionBordXm == 4))) ||
            ((!strcmp(TypeBord,"Ym")) && ((ConditionBordXm == 2)||(ConditionBordXm == 4))))
    {
		#pragma omp parallel for private(Legerete) schedule(guided)
        for(i = 1 ; i < W ; i++)
        {
            Legerete = (TableLegerete[Carte_Tampon[0]]);
            Vj_I[i*(J+1)+0] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                                  (1-DeltaT*AttenuationPML[2*i])*Vj_I[i*(J+1)+0]+ DeltaTsurH*Legerete*SignI*(
                                      (Tij_I[(i)*(J+1)+0] + Tij_J[(i)*(J+1)+0]) - (Tij_I[(i-1)*(J+1)+0] + Tij_J[(i-1)*(J+1)+0]))) ;
            Vj_J[i*(J+1)+0] += DeltaTsurH*Legerete*(
                                   2*(Tjj_I[i*J+0]	+ Tjj_J[i*J+0]) - 0) ;
        }
        /*Coin du bord PML*/
        Legerete = (TableLegerete[Carte_Tampon[0]]) ;
        Vj_I[0*(J+1)+0] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                              (1-DeltaT*AttenuationPML[2*0])*Vj_I[0*(J+1)+0]+ DeltaTsurH*Legerete*SignI*(
                                  (Tij_I[(0)*(J+1)+0] + Tij_J[(0)*(J+1)+0]) - (Tij_BufferBordBlocI[0]))) ;
        Vj_J[0*(J+1)+0] += DeltaTsurH*Legerete*(
                               2*(Tjj_I[0*J+0]	+ Tjj_J[0*J+0]) - 0) ;
    }
    /**************************************************************************/
    /*********************** Attenuating Factor (QModel) **********************/
    /**************************************************************************/
    if(TypeAbsorption == 2)
    {
		#pragma omp parallel for private(j,FacteurQModel) schedule(guided)
        for(i = 0 ; i < (W) ; i++)
        {
            for(j = 0 ; j < J ; j++)
            {
                FacteurQModel = TableQModel[Carte_Tampon[j]] ;
                Vi_I[i*(J)+j] *= FacteurQModel ;
                Vi_J[i*(J)+j] *= FacteurQModel ;
            }
        }
		#pragma omp parallel for private(j,FacteurQModel) schedule(guided)
        for(i = 0 ; i < (W+1) ; i++)
        {
            for(j = 1 ; j < J ; j++)
            {
                FacteurQModel = sqrt(TableQModel[Carte_Tampon[j-1]]*TableQModel[Carte_Tampon[j]]) ;
                Vj_I[i*(J+1)+j] *= FacteurQModel ;
                Vj_J[i*(J+1)+j] *= FacteurQModel ;
            }
            FacteurQModel = TableQModel[Carte_Tampon[j]] ;
            Vj_I[i*(J+1)+J] *= FacteurQModel ;
            FacteurQModel = TableQModel[Carte_Tampon[0]] ;
            Vj_I[i*(J+1)+0] *= FacteurQModel ;
        }
    }
}
