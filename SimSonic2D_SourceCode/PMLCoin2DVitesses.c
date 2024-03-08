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

/** \file PMLCoin2DVitesses.c
    \brief Code source des fonctions de calcul des PML pour les vitesses dans les coins de la carte.
    \author Emmanuel Bossy
    \version 2D_2011_10_28_omp
    \date 28 octobre 2011 */

void CalculPMLCoinVitesses(PREC_GRILLE *Vi_I,PREC_GRILLE *Vi_J,
                           PREC_GRILLE *Vj_I,PREC_GRILLE *Vj_J, PREC_GRILLE *Tii_I,PREC_GRILLE *Tii_J,
                           PREC_GRILLE *Tjj_I,PREC_GRILLE *Tjj_J, PREC_GRILLE *Tij_I,PREC_GRILLE *Tij_J,
                           const char *TypeCoin)
{
    int i,j ;
    double Legerete ;
    double FacteurQModel ;
    double SignI, SignJ ;

    if (!strcmp(TypeCoin,"XpYp") || !strcmp(TypeCoin,"XpYm"))
        SignI = +1.0 ;
    else
        SignI = -1.0 ;

    if (!strcmp(TypeCoin,"XpYp") || !strcmp(TypeCoin,"XmYp"))
        SignJ = +1.0 ;
    else
        SignJ = -1.0 ;

    Legerete = TableLegerete[Carte_Tampon[0]] ;
    /*Bloc...*/
	#pragma omp parallel for private(j) schedule(guided)
    for(i = 0 ; i < W ; i++)
    {
        for(j = 1 ; j < W ; j++)
        {
            Vi_I[i*(W+1)+j] = 1./(1+DeltaT*AttenuationPML[1+2*i])*(
                                  (1-DeltaT*AttenuationPML[1+2*i])*Vi_I[i*(W+1)+j]+ DeltaTsurH*Legerete*SignI*(
                                      (Tii_I[(i+1)*(W+1)+j] + Tii_J[(i+1)*(W+1)+j]) - (Tii_I[(i)*(W+1)+j] + Tii_J[(i)*(W+1)+j]))) ;
            Vi_J[i*(W+1)+j] = 1./(1+DeltaT*AttenuationPML[2*j])*(
                                  (1-DeltaT*AttenuationPML[2*j])*Vi_J[i*(W+1)+j]+ DeltaTsurH*Legerete*SignJ*(
                                      (Tij_I[i*(W)+j]	+ Tij_J[i*(W)+j]) - (Tij_I[i*(W)+j-1] +	Tij_J[i*(W)+j-1]))) ;
        }
    }
	#pragma omp parallel for private(j) schedule(guided)
    for(i = 1 ; i < W ; i++)
    {
        for(j = 0 ; j < W ; j++)
        {
            Vj_I[i*(W)+j] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                                (1-DeltaT*AttenuationPML[2*i])*Vj_I[i*(W)+j]+ DeltaTsurH*Legerete*SignI*(
                                    (Tij_I[(i)*(W)+j] + Tij_J[(i)*(W)+j]) -	(Tij_I[(i-1)*(W)+j] + Tij_J[(i-1)*(W)+j]))) ;
            Vj_J[i*(W)+j] = 1./(1+DeltaT*AttenuationPML[1+2*j])*(
                                (1-DeltaT*AttenuationPML[1+2*j])*Vj_J[i*(W)+j]+ DeltaTsurH*Legerete*SignJ*(
                                    (Tjj_I[i*(W+1)+j+1] + Tjj_J[i*(W+1)+j+1]) - (Tjj_I[i*(W+1)+j] +	Tjj_J[i*(W+1)+j]))) ;
        }
    }
    /*Frontières*/
	#pragma omp parallel for schedule(guided)
    for(i = 0 ; i < W ; i++)
    {
        Vi_I[i*(W+1)+0] = 1./(1+DeltaT*AttenuationPML[1+2*i])*(
                              (1-DeltaT*AttenuationPML[1+2*i])*Vi_I[i*(W+1)+0]+ DeltaTsurH*Legerete*SignI*(
                                  (Tii_I[(i+1)*(W+1)+0] +	Tii_J[(i+1)*(W+1)+0]) -	(Tii_I[(i)*(W+1)+0] + Tii_J[(i)*(W+1)+0]))) ;
        Vi_J[i*(W+1)+0] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                              (1-DeltaT*AttenuationPML[2*0])*Vi_J[i*(W+1)+0]+ DeltaTsurH*Legerete*SignJ*(
                                  (Tij_I[i*(W)+0]	+ Tij_J[i*(W)+0]) - (Tij_BufferCoinBordI[i]))) ;
    }
	#pragma omp parallel for schedule(guided)
    for(j = 0 ; j < W ; j++)
    {
        Vj_I[0*(W)+j] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                            (1-DeltaT*AttenuationPML[2*0])*Vj_I[0*(W)+j]+ DeltaTsurH*Legerete*SignI*(
                                (Tij_I[(0)*(W)+j] + Tij_J[(0)*(W)+j]) - (Tij_BufferCoinBordJ[j]))) ;
        Vj_J[0*(W)+j] = 1./(1+DeltaT*AttenuationPML[1+2*j])*(
                            (1-DeltaT*AttenuationPML[1+2*j])*Vj_J[0*(W)+j]+ DeltaTsurH*Legerete*SignJ*(
                                (Tjj_I[0*(W+1)+j+1] + Tjj_J[0*(W+1)+j+1]) - (Tjj_I[0*(W+1)+j] +	Tjj_J[0*(W+1)+j]))) ;
    }
    /**************************************************************************/
    /************************Attenuating Factor (QModel)***********************/
    /**************************************************************************/
    if(TypeAbsorption == 2)
    {
        FacteurQModel = TableQModel[Carte_Tampon[0]];
		#pragma omp parallel for private(j) schedule(guided)
        for(i = 0 ; i < W ; i++)
        {
            for(j = 0 ; j < (W+1) ; j++)
            {
                Vi_I[i*(W+1)+j] *= FacteurQModel ;
                Vi_J[i*(W+1)+j] *= FacteurQModel ;
            }
        }
		#pragma omp parallel for private(j) schedule(guided)
        for(i = 0 ; i < (W+1) ; i++)
        {
            for(j = 0 ; j < W ; j++)
            {
                Vj_I[i*(W)+j] *= FacteurQModel ;
                Vj_J[i*(W)+j] *= FacteurQModel ;
            }
        }
    }
}
