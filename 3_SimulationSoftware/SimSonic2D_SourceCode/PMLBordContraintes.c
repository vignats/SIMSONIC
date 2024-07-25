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


/** \file PMLBordContraintes.c
    \brief Code source des fonctions de calcul des PML pour les contraintes sur les bords de la carte.
    \author Emmanuel Bossy
    \version 2D_2011_10_28_omp
    \date 28 octobre 2011 */

void CalculPMLBordContraintes(PREC_GRILLE *Vi_I,PREC_GRILLE *Vi_J,
                              PREC_GRILLE *Vj_I,PREC_GRILLE *Vj_J, PREC_GRILLE *Tii_I,PREC_GRILLE *Tii_J,
                              PREC_GRILLE *Tjj_I,PREC_GRILLE *Tjj_J, PREC_GRILLE *Tij_I,PREC_GRILLE *Tij_J,
                              int J, const char *TypeBord)
{
    int i,j ;
    double SignI ;
    double Cijij ;
    double	FacteurQModel ;
    if(!strcmp(TypeBord,"Xp") || !strcmp(TypeBord,"Yp"))
        SignI = +1.0 ;
    else
        SignI = -1.0 ;
    /*************************************************************************/
    /***************************  Bloc du Bord PML ***************************/
    /*************************************************************************/


	#pragma omp parallel for private(j) schedule(guided)
    for(i = 1 ; i < (W) ; i++)
    {
        for(j = 0 ; j < J ; j++)
        {
            Tii_I[i*J+j] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                               (1-DeltaT*AttenuationPML[2*i])*Tii_I[i*J+j] +
                               DeltaTsurH*TableCii[Carte_Tampon[j]]*SignI*(
                                   (Vi_I[(i)*J+j] + Vi_J[(i)*J+j]) - (Vi_I[(i-1)*J+j] + Vi_J[(i-1)*J+j]))) ;
            Tii_J[i*J+j] += DeltaTsurH*TableCij[Carte_Tampon[j]]*
                            ((Vj_I[i*(J+1)+j+1] + Vj_J[i*(J+1)+j+1]) - (Vj_I[i*(J+1)+j] + Vj_J[i*(J+1)+j])) ;
            Tjj_I[i*J+j] = 1./(1+DeltaT*AttenuationPML[2*i])*(
                               (1-DeltaT*AttenuationPML[2*i])*Tjj_I[i*J+j] +
                               DeltaTsurH*TableCij[Carte_Tampon[j]]*SignI*(
                                   (Vi_I[(i)*J+j] + Vi_J[(i)*J+j]) - (Vi_I[(i-1)*J+j] + Vi_J[(i-1)*J+j]))) ;
            Tjj_J[i*J+j] += DeltaTsurH*TableCjj[Carte_Tampon[j]]*
                            ((Vj_I[i*(J+1)+j+1] + Vj_J[i*(J+1)+j+1]) - (Vj_I[i*(J+1)+j] + Vj_J[i*(J+1)+j])) ;
        }
    }



    #pragma omp parallel for private(j,Cijij) schedule(guided)
	for(i=0;i<(W);i++)
    {
        for(j=1;j<J;j++)
        {
            if(TableCijij[Carte_Tampon[j-1]]*TableCijij[Carte_Tampon[j]])
                Cijij=2./(1./TableCijij[Carte_Tampon[j-1]]+1./TableCijij[Carte_Tampon[j]]) ;
            else
                Cijij=0.0;
            Tij_I[i*(J+1)+j]=1./(1+DeltaT*AttenuationPML[1+2*i])*(
                                   (1-DeltaT*AttenuationPML[1+2*i])*Tij_I[i*(J+1)+j] +  DeltaTsurH*SignI*Cijij*(
                                       (Vj_I[(i+1)*(J+1)+j] + Vj_J[(i+1)*(J+1)+j]) - (Vj_I[(i)*(J+1)+j] + Vj_J[(i)*(J+1)+j])));
            Tij_J[i*(J+1)+j] += DeltaTsurH*Cijij*(
                                    (Vi_I[i*J+j] + Vi_J[i*J+j]) - (Vi_I[i*J+j-1] + Vi_J[i*J+j-1]));
        }
    }
    /*************************************************************************/
    /*************************** Bord sur Bloc du Bord PML ...****************/
    /*************************************************************************/



	#pragma omp parallel for schedule(guided)
    for(j = 0 ; j < J ; j++)
    {
        Tii_I[0*J+j] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                           (1-DeltaT*AttenuationPML[2*0])*Tii_I[0*J+j] +
                           DeltaTsurH*TableCii[Carte_Tampon[j]]*SignI*(
                               (Vi_I[(0)*J+j] + Vi_J[(0)*J+j]) - (Vi_BufferBordBlocI[j]))) ;
        Tii_J[0*J+j] += DeltaTsurH*TableCij[Carte_Tampon[j]]*
                        ((Vj_I[0*(J+1)+j+1] + Vj_J[0*(J+1)+j+1]) - (Vj_I[0*(J+1)+j] + Vj_J[0*(J+1)+j])) ;
        Tjj_I[0*J+j] = 1./(1+DeltaT*AttenuationPML[2*0])*(
                           (1-DeltaT*AttenuationPML[2*0])*Tjj_I[0*J+j] +
                           DeltaTsurH*TableCij[Carte_Tampon[j]]*SignI*(
                               (Vi_I[(0)*J+j] + Vi_J[(0)*J+j]) - (Vi_BufferBordBlocI[j]))) ;
        Tjj_J[0*J+j] += DeltaTsurH*TableCjj[Carte_Tampon[j]]*
                        ((Vj_I[0*(J+1)+j+1] + Vj_J[0*(J+1)+j+1]) - (Vj_I[0*(J+1)+j] + Vj_J[0*(J+1)+j])) ;
    }
    /*************************************************************************/
    /*************************** Bord sur Coin du Bord PML ...****************/
    /*************************************************************************/
    if(((!strcmp(TypeBord,"Xp")) && (ConditionBordYm == 0)) ||
            ((!strcmp(TypeBord,"Xm")) && (ConditionBordYm == 0)) ||
            ((!strcmp(TypeBord,"Yp")) && (ConditionBordXm == 0)) ||
            ((!strcmp(TypeBord,"Ym")) && (ConditionBordXm == 0)))
    {


		#pragma omp parallel for private(Cijij) schedule(guided)
        for(i = 0 ; i < (W) ; i++)
        {
            Cijij = TableCijij[Carte_Tampon[0]] ;
            Tij_I[i*(J+1)+0] = 1./(1+DeltaT*AttenuationPML[1+2*i])*(
                                   (1-DeltaT*AttenuationPML[1+2*i])*Tij_I[i*(J+1)+0] +  DeltaTsurH*SignI*Cijij*(
                                       (Vj_I[(i+1)*(J+1)+0] + Vj_J[(i+1)*(J+1)+0]) - (Vj_I[(i)*(J+1)+0] + Vj_J[(i)*(J+1)+0]))) ;
            Tij_J[i*(J+1)+0] += DeltaTsurH*Cijij*(
                                    (Vi_I[i*J+0] + Vi_J[i*J+0]) - (Vi_BufferBordCoinJm[i])) ;
        }
    }
    if(((!strcmp(TypeBord,"Xp")) && (ConditionBordYp == 0)) ||
            ((!strcmp(TypeBord,"Xm")) && (ConditionBordYp == 0)) ||
            ((!strcmp(TypeBord,"Yp")) && (ConditionBordXp == 0)) ||
            ((!strcmp(TypeBord,"Ym")) && (ConditionBordXp == 0)))
    {
		#pragma omp parallel for private(Cijij) schedule(guided)
        for(i = 0 ; i < (W) ; i++)
        {
            Cijij = TableCijij[Carte_Tampon[J-1]] ;
            Tij_I[i*(J+1)+J] = 1./(1+DeltaT*AttenuationPML[1+2*i])*(
                                   (1-DeltaT*AttenuationPML[1+2*i])*Tij_I[i*(J+1)+J] +  DeltaTsurH*SignI*Cijij*(
                                       (Vj_I[(i+1)*(J+1)+J] + Vj_J[(i+1)*(J+1)+J]) -(Vj_I[(i)*(J+1)+J] + Vj_J[(i)*(J+1)+J]))) ;
            Tij_J[i*(J+1)+J] += DeltaTsurH*Cijij*(
                                    (Vi_BufferBordCoinJp[i]) - (Vi_I[i*J+J-1] + Vi_J[i*J+J-1])) ;
        }
    }

    /***************** Frontières Coin Rigide...*********************************/
    if(((!strcmp(TypeBord,"Xp")) && ((ConditionBordYm == 3)||(ConditionBordYm == 4))) ||
            ((!strcmp(TypeBord,"Xm")) && ((ConditionBordYm == 3)||(ConditionBordYm == 4))) ||
            ((!strcmp(TypeBord,"Yp")) && ((ConditionBordXm == 3)||(ConditionBordXm == 4))) ||
            ((!strcmp(TypeBord,"Ym")) && ((ConditionBordXm == 3)||(ConditionBordXm == 4))))
    {
		#pragma omp parallel for private(Cijij) schedule(guided)
        for(i = 0 ; i < (W) ; i++)
        {
            Cijij = TableCijij[Carte_Tampon[0]] ;
            Tij_I[i*(J+1)+0] = 1./(1+DeltaT*AttenuationPML[1+2*i])*(
                                   (1-DeltaT*AttenuationPML[1+2*i])*Tij_I[i*(J+1)+0] +  DeltaTsurH*SignI*Cijij*(
                                       (Vj_I[(i+1)*(J+1)+0] + Vj_J[(i+1)*(J+1)+0]) - (Vj_I[(i)*(J+1)+0] + Vj_J[(i)*(J+1)+0]))) ;
            Tij_J[i*(J+1)+0] += DeltaTsurH*Cijij*(
                                    2*(Vi_I[i*J+0] + Vi_J[i*J+0]) -	0) ;
        }
    }
    if(((!strcmp(TypeBord,"Xp")) && ((ConditionBordYp == 3)||(ConditionBordYp == 4))) ||
            ((!strcmp(TypeBord,"Xm")) && ((ConditionBordYp == 3)||(ConditionBordYp == 4))) ||
            ((!strcmp(TypeBord,"Yp")) && ((ConditionBordXp == 3)||(ConditionBordXp == 4))) ||
            ((!strcmp(TypeBord,"Ym")) && ((ConditionBordXp == 3)||(ConditionBordXp == 4))))
    {
		#pragma omp parallel for private(Cijij) schedule(guided)
        for(i = 0 ; i < (W) ; i++)
        {
            Cijij = TableCijij[Carte_Tampon[J-1]] ;
            Tij_I[i*(J+1)+J] = 1./(1+DeltaT*AttenuationPML[1+2*i])*(
                                   (1-DeltaT*AttenuationPML[1+2*i])*Tij_I[i*(J+1)+J] +  DeltaTsurH*SignI*Cijij*(
                                       (Vj_I[(i+1)*(J+1)+J] + Vj_J[(i+1)*(J+1)+J]) - (Vj_I[(i)*(J+1)+J] + Vj_J[(i)*(J+1)+J]))) ;
            Tij_J[i*(J+1)+J] += DeltaTsurH*Cijij*(
                                    0 - 2*(Vi_I[i*J+J-1] + Vi_J[i*J+J-1])) ;
        }
    }
    /**************************************************************************/
    /*********************** Attenuating Factor (QModel) **********************/
    /**************************************************************************/
    if(TypeAbsorption == 2)
    {
		#pragma omp parallel for private(j,FacteurQModel) schedule(guided)
        for(i = 0 ; i < (W+1) ; i++)
        {
            for(j = 0 ; j < J ; j++)
            {
                FacteurQModel = TableQModel[Carte_Tampon[j]] ;
                Tii_I[i*J+j] *= FacteurQModel ;
                Tii_J[i*J+j] *= FacteurQModel ;
                Tjj_I[i*J+j] *= FacteurQModel ;
                Tjj_J[i*J+j] *= FacteurQModel ;
            }
        }
		#pragma omp parallel for private(j,FacteurQModel) schedule(guided)
        for(i = 0 ; i < (W) ; i++)
        {
            for(j = 1 ; j < J ; j++)
            {
                FacteurQModel = sqrt(TableQModel[Carte_Tampon[j-1]]*TableQModel[Carte_Tampon[j]]) ;
                Tij_I[i*(J+1)+j] *= FacteurQModel ;
                Tij_J[i*(J+1)+j] *= FacteurQModel ;
            }
            FacteurQModel = TableQModel[Carte_Tampon[J-1]] ;
            Tij_J[i*(J+1)+J] *= FacteurQModel ;
            FacteurQModel = TableQModel[Carte_Tampon[0]] ;
            Tij_J[i*(J+1)+0] *= FacteurQModel ;
        }
    }
}
