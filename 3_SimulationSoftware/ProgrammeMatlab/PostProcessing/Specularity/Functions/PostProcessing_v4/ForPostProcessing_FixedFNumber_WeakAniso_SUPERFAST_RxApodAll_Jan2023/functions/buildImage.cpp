#include "buildImage.hpp"
using namespace std;

void*
buildImage(const ExpSetup_st *Setup,const Input_st *Data, Output_st *Image)
{
    int i, j, iTx, Tx, r, Index, OffsetT, OffsetSig, OffsetR;
    int min_ap1_idx, ap1_count, min_ap2_idx, ap2_count, min_ap3_idx, ap3_count, min_ap4_idx, ap4_count, previous_element;
    double TotalTime, Weight, TransmitTime, compare_ap; 

    omp_set_num_threads(Setup->NCore);
    #pragma omp parallel for schedule(static) \
    shared(Setup, Data, Image) \
    private( i,j, iTx, Tx, OffsetT, TransmitTime, compare_ap, \
                    r, OffsetR, TotalTime, Index, \
                    OffsetSig, Weight, min_ap1_idx, ap1_count , \
                    min_ap2_idx, ap2_count, min_ap3_idx, ap3_count, min_ap4_idx, ap4_count , previous_element )
    for (i = 0; i < Data->NX; ++i)
	{
        double* pWindow = new double[Data->NR]; // modified by GR, Oct 2022
        
        for (j = Data->ZMin; j < Data->ZMax; ++j)
        {
            // Receive apodization, modified by GR, Oct 2022
            min_ap1_idx = 0;
            ap1_count = 0;
            min_ap2_idx = 0;
            ap2_count = 0;
            min_ap3_idx = 0;
            ap3_count = 0;
            min_ap4_idx = 0;
            ap4_count = 0;
            previous_element = 0;
            for (r = 0; r < Data->NR; ++r)
            {
                OffsetR = i * Data->NR * Data->NZ + j * Data->NR + r;

                if (fabs(Image->pAngleR[OffsetR]) < Setup->RxHalfOpenAngle & Image->pTimeR[OffsetR] < 1.) // travel_time=1 if no physical ray
                {
                    if (previous_element == 0) // start new subaperture
                    {
                    if (ap1_count == 0) min_ap1_idx = r; // first element of first receive subaperture
                    else if (ap2_count == 0) min_ap2_idx = r; // first element of second receive subaperture
                    else if (ap3_count == 0) min_ap3_idx = r; // first element of third receive subaperture
                    else min_ap4_idx = r; // first element of fourth receive subaperture
                    }
                    
                    if (min_ap2_idx == 0) ap1_count++; // increment number of elements in first receive subaperture
                    else if (min_ap3_idx == 0) ap2_count++; // increment number of elements in second receive subaperture
                    else if (min_ap4_idx == 0) ap3_count++; // increment number of elements in third receive subaperture
                    else ap4_count++; // increment number of elements in fourth receive subaperture
                    
                    previous_element = 1;
                }
                else
                {
                    previous_element = 0;
                }
            }
            
            
            if ( ap2_count > 0 ) // if two subapertures, select largest and nearest subaperture to the pixel
            {
            compare_ap = ap2_count/ap1_count * fabs(Data->pX[i]-Data->pXR[int(min_ap1_idx+ap1_count/2.)])/fabs(Data->pX[i]-Data->pXR[int(min_ap2_idx+ap2_count/2.)]);
                if (compare_ap > 1.)
                {
                    ap1_count = ap2_count;
                    min_ap1_idx = min_ap2_idx;
                }
            }
            
            if ( ap3_count > 0 ) // if three subapertures, select largest and nearest subaperture to the pixel
            {
            compare_ap = ap3_count/ap1_count * fabs(Data->pX[i]-Data->pXR[int(min_ap1_idx+ap1_count/2.)])/fabs(Data->pX[i]-Data->pXR[int(min_ap3_idx+ap3_count/2.)]);
                if (compare_ap > 1.)
                {
                    ap1_count = ap3_count;
                    min_ap1_idx = min_ap3_idx;
                }
            }
            
            if ( ap4_count > 0 ) // if four subapertures, select largest and nearest subaperture to the pixel
            {
            compare_ap = ap4_count/ap1_count * fabs(Data->pX[i]-Data->pXR[int(min_ap1_idx+ap1_count/2.)])/fabs(Data->pX[i]-Data->pXR[int(min_ap4_idx+ap4_count/2.)]);
                if (compare_ap > 1.)
                {
                    ap1_count = ap4_count;
                    min_ap1_idx = min_ap4_idx;
                }
            }

            
            if ( (min_ap1_idx+ap1_count) > (Data->NR-1) ) ap1_count = (Data->NR-1) - min_ap1_idx; // safety in case more than 4 subapertures     

            
            
            if (ap1_count > 0)
            {
                
            setArrayToZero(pWindow, ap1_count);
            getWindowFunction(pWindow, ap1_count, Setup->SubAperApodis);
                
                for (r = min_ap1_idx; r < (min_ap1_idx+ap1_count); ++r) // modified by GR, Nov 2022
                {
                OffsetR = i*Data->NR*Data->NZ + j*Data->NR + r;
                if (fabs(Image->pAngleR[OffsetR]) < Setup->RxHalfOpenAngle)
                { // if angle at receive element lies within acceptance angle
                    for (iTx = 0; iTx < Data->NTX; ++iTx)  // loop on sources
                    {
                        Tx 	= Data->pTxToUse[iTx] - 1;
                        OffsetT = i*Data->NS*Data->NZ + j*Data->NS + Tx;
                        TransmitTime = Image->pTimeT[OffsetT]  - Data->pAddDelay[Tx];
                        
                        if (fabs(Image->pAngleT[OffsetT]) < Setup->TxHalfOpenAngle & TransmitTime < 1.) // travel_time=1 if no physical ray
                        { // if angle at transmit element lies within acceptance angle
                            TotalTime = (Image->pTimeR[OffsetR] + TransmitTime) * Setup->SamplingFreq;
                            Index     = int(TotalTime);

                            if (Index >= 1 && Index < Data->TimeSig -1)
                            { // Delay and Sum, linear interpolation
                                OffsetSig = Tx*Data->TimeSig*Data->NR + r*Data->TimeSig + Index;
                                Weight    = 1. - (TotalTime - floor(TotalTime));
                                Image->pImI[i][j]     += (Data->pSignalI[OffsetSig] * Weight + Data->pSignalI[OffsetSig + 1] * (1. - Weight)) * pWindow[r - min_ap1_idx];
                                Image->pImQ[i][j]     += (Data->pSignalQ[OffsetSig] * Weight + Data->pSignalQ[OffsetSig + 1] * (1. - Weight)) * pWindow[r - min_ap1_idx];

                            }
                        }
                    } // end loop on sources
                }
            } // end loop on receivers
        }
        } // end loop on z (image depth) dimension
        delete[] pWindow;  // modified by GR, Sept 2022        
    } // end loop on x (image width) dimension
    return NULL;
}

void*
buildImageForVFI(const ExpSetup_st *Setup, const Input_st *Data, Output_st *Image)
{
    int i, j, iTx, Tx, iRx, r, Index, OffsetT, OffsetSig, OffsetR, OffsetBeam;
    double sumApI, sumApQ, TotalTime, Weight, TempI, TempQ, TransmitTime;
    int minR, maxR, nElemInSubR;

    omp_set_num_threads(Setup->NCore);
    #pragma omp parallel for schedule(static) \
    shared(Setup, Data, Image) \
    private(i, j, iTx, Tx, OffsetT, TransmitTime, iRx, sumApI, sumApQ, r, OffsetR, \
            TotalTime, Index, OffsetSig, Weight, TempI, TempQ, OffsetBeam, minR, maxR, nElemInSubR )
    for (i = 0; i < Data->NX; ++i)
	{
        double* pWindow = new double [Data->NR];
        for (j = Data->ZMin; j < Data->ZMax; ++j)
        {
            for (iTx = 0; iTx < Data->NTX; ++iTx)
            { // Nb of Transmission.
                Tx 	= int(Data->pTxToUse[iTx]) - 1;
                OffsetT = i*Data->NS*Data->NZ + j*Data->NS + Tx;
                
                if (fabs(Image->pAngleT[OffsetT]) < Setup->TxHalfOpenAngle & Image->pTimeT[OffsetT] < 1.) // travel_time=1 if no physical ray
                { // if angle at transmit element lies within acceptance angle

                TransmitTime = Image->pTimeT[OffsetT] - Data->pAddDelay[Tx];

                for (iRx = 0; iRx < Data->NLA; ++iRx)
                { // Nb of listening angles
                    sumApI = 0.; sumApQ = 0.; minR = 0.; maxR = 0.; nElemInSubR = 0;
                    // Search for the active listening subaperture for a given listening angle anf F-number.
                    nElemInSubR = getActiveSubAperture(i, j, iRx, Image, Data, &minR, &maxR, Setup->Fnumber, Setup->MaxDevListenAngle);
                    // Initialize Window at each iRx
                    //setArrayToZero(pWindow,Data->NR);
                    setArrayToZero(pWindow,nElemInSubR); // modified by GR, October 2022
                    // If we have listening elements (size of subaperture!=0) we build.
                    if (nElemInSubR > 1) // modified by GR, January 2023
                    {
                        // Apply apodization on listening SubAperture.
                        getWindowFunction(pWindow, nElemInSubR, Setup->SubAperApodis);
                        
                        for (r = minR; r <= maxR; r++)
                        { // receiver index within Rx aperture
                            OffsetR   = i*Data->NR*Data->NZ + j*Data->NR + r;
                            TotalTime = (Image->pTimeR[OffsetR] + TransmitTime) * Setup->SamplingFreq;
                            Index     = int(TotalTime);
                            
                            if (Index >= 1 && Index < Data->TimeSig -1)
                            { // Delay and Sum, linear interpolation
                                OffsetSig = Tx*Data->TimeSig*Data->NR + r*Data->TimeSig + Index;
                                Weight    = 1. - (TotalTime - floor(TotalTime));
                                TempI     = Data->pSignalI[OffsetSig] * Weight + Data->pSignalI[OffsetSig + 1] * (1. - Weight);
                                TempQ     = Data->pSignalQ[OffsetSig] * Weight + Data->pSignalQ[OffsetSig + 1] * (1. - Weight);
                                sumApI    += TempI*pWindow[r-minR];
                                sumApQ    += TempQ*pWindow[r-minR];
                            }
                        }
                        OffsetBeam = Tx*Data->NX*Data->NZ*Data->NLA +
                                        iRx*Data->NX*Data->NZ +
                                        i*Data->NZ +
                                        j;

                    Image->pBeamQ[OffsetBeam] += sumApQ;
                    Image->pBeamI[OffsetBeam] += sumApI;
                    }
                } // end loop on receive angles
            }
            } // end loop on transmissions
        } // end loop on z (image depth) dimension
        delete [] pWindow;
    } // end loop on x (image width) dimension
    return NULL;
}

int
getActiveSubAperture(const int IdX, const int IdZ, const int IdL, const Output_st* pImage, const Input_st* pData,
                     int *rMinR, int *rMaxR, const double FNumber, const double MaxDeviationFromListenAngle)
{
    double diffAngle, minDiffAngle = 1.0; // modified by GR, Nov 2022, 0.5236; // Difference between LisenAngle/AngleAtPix, Initial value for minimum seeking. (30deg)
    int OffsetR, IdxClosest0;
    double IdxClosest                 = 0.0;

    // We first get the element which is aligned with the given listenAngle from a given Pixel(x,z).
    for (int r = 0; r < pData->NR; r++)
    {// Looping over elements to find r which observe a pixel with the closest value to ListenAngle.

        OffsetR = IdX*pData->NR*pData->NZ + IdZ*pData->NR + r; // added by GR, Nov 2022
        for (double half_el = 0; half_el < 1; half_el = half_el + 0.5)
        {
            if (half_el == 0) // central ray hits a physical array element
            {
            // Difference between the angle of refraction of the wave and the ListenAngle.
            diffAngle = fabs(pData->pListenAngle[IdL] - pImage->pAngleAtPixR[OffsetR]); // modified by GR, Nov 2022
                if ( diffAngle < minDiffAngle & pImage->pTimeR[OffsetR] < 1.) // travel_time=1 if no physical ray
                {
                        minDiffAngle = diffAngle; // keep minimum difference.
                        IdxClosest   = r;         // IdxCloset = index of r that minimzes.
                }
            }
            else if (r < pData->NR-1) // central ray hits the center between two physical array elements
            {
            // Difference between the angle of refraction of the wave and the ListenAngle.
            diffAngle = fabs(pData->pListenAngle[IdL] - (pImage->pAngleAtPixR[OffsetR]+pImage->pAngleAtPixR[OffsetR+1])/2); // added by GR, Jan 2023
                if ( diffAngle < minDiffAngle & pImage->pTimeR[OffsetR] < 1.) // travel_time=1 if no physical ray
                {
                        minDiffAngle = diffAngle; // keep minimum difference.
                        IdxClosest   = r + 0.5;         // IdxCloset = index of r that minimzes.
                }
            }
        }
    }
    
    // If there are no element listening with that angle (Or too far from that listening angle, i.e.  > tresholdAngle, we dont reconstruct image.
    if (fabs(minDiffAngle) > MaxDeviationFromListenAngle) return 0; // Not very useful, but for safety.

    // Else IdxClosest is the center of an active listening subaperture.
    // Using the F-number definition and its value imposed by user, we get the width of this subaperture centered at IdxClosest.
    // Approximation: The subaperture is symmetrical around IdxClosest.
    else {
        // We keep angles at Pixel when the center is identified.
        if (IdxClosest-floor(IdxClosest) == 0) // central ray hits a physical array element
        {
            IdxClosest0 = int(IdxClosest);
            pImage->pListenedAngleAtPixR[IdX*pData->NLA*pData->NZ + IdZ*pData->NLA + IdL] =
                    pImage->pAngleAtPixR[IdX*pData->NR*pData->NZ + IdZ*pData->NR + IdxClosest0];

            // Euclidian distance between center of sub and a given pixel.
            double FocalDistance       = distance(pData->pXR[IdxClosest0], 0, pData->pX[IdX], pData->pZ[IdZ]);
            double ActiveApertureWidth = FocalDistance/FNumber;

            // Get the 1st element within Activeaperturewidth left to Idxclosest;
            for (int i = 0; i <= IdxClosest0; i++)// modified by GR, Nov 2022, <= instead of <
            {
                if (fabs(pData->pXR[i]-pData->pXR[IdxClosest0]) >= ActiveApertureWidth/2) *rMinR = i; // modified by GR, Nov 2022, fabs instead of sqrt(sqr(
            }

            // Get the 1st element within Activeaperturewidth right to Idxclosest;
            for (int i = IdxClosest0; i < pData->NR; i++)
            {
                if (fabs(pData->pXR[i]-pData->pXR[IdxClosest0]) <= ActiveApertureWidth/2) *rMaxR = i; // modified by GR, Nov 2022, fabs instead of sqrt(sqr(
            }

            // Check if the subaperture is symmetrical. If yes we move.
            if (*rMaxR-IdxClosest0==IdxClosest0-*rMinR) return abs(*rMaxR - *rMinR) + 1;

            // if not we check who is the longest and scale down to the other side value.
            if (*rMaxR-IdxClosest0 > IdxClosest0-*rMinR) *rMaxR = 2*IdxClosest0 - *rMinR; // because we want: *rMaxR - IdxClosest = IdxClosest - *rMinR
            else *rMinR = 2*IdxClosest0 - *rMaxR;

            // Number of element within our Sub.
            return abs(*rMaxR - *rMinR) + 1;
        }
        else  // central ray hits the center between two physical array elements
        {
            pImage->pListenedAngleAtPixR[IdX*pData->NLA*pData->NZ + IdZ*pData->NLA + IdL] =
                    (pImage->pAngleAtPixR[IdX*pData->NR*pData->NZ + IdZ*pData->NR + int(IdxClosest)]+
                    pImage->pAngleAtPixR[IdX*pData->NR*pData->NZ + IdZ*pData->NR + int(IdxClosest)+1])/2;

            // Euclidian distance between center of sub and a given pixel.
            double XCenterAperture     = (pData->pXR[int(IdxClosest)]+pData->pXR[int(IdxClosest)+1])/2;
            double FocalDistance       = distance(XCenterAperture, 0, pData->pX[IdX], pData->pZ[IdZ]);
            double ActiveApertureWidth = FocalDistance/FNumber;

            // Get the 1st element within Activeaperturewidth left to Idxclosest;
            for (int i = 0; i <= int(floor(IdxClosest)); i++)// modified by GR, Nov 2022, <= instead of <
            {
                if (fabs(pData->pXR[i]-XCenterAperture) >= ActiveApertureWidth/2) *rMinR = i; // modified by GR, Nov 2022, fabs instead of sqrt(sqr(
            }

            // Get the 1st element within Activeaperturewidth right to Idxclosest;
            for (int i = int(ceil(IdxClosest)); i < pData->NR; i++)
            {
                if (fabs(pData->pXR[i]-XCenterAperture) <= ActiveApertureWidth/2) *rMaxR = i; // modified by GR, Nov 2022, fabs instead of sqrt(sqr(
            }

            // Check if the subaperture is symmetrical. If yes we move.
            if (double(*rMaxR)-IdxClosest==IdxClosest-double(*rMinR)) return abs(*rMaxR - *rMinR) + 1;

            // if not we check who is the longest and scale down to the other side value.
            if (double(*rMaxR)-IdxClosest > IdxClosest-double(*rMinR)) *rMaxR = int(2*IdxClosest) - *rMinR; // because we want: *rMaxR - IdxClosest = IdxClosest - *rMinR
            else *rMinR = int(2*IdxClosest) - *rMaxR;

            // Number of element within our Sub.
            return abs(*rMaxR - *rMinR) + 1;
        }
    }
}

double
getFnumber(const double ElemXMin, const double ElemXMax, const double PixelX, const double PixelZ)
{
    return distance((ElemXMax+ElemXMin)/2, 0, PixelX, PixelZ) / (ElemXMax-ElemXMin);
}

void*
getWindowFunction(double* pWindow,  const int N, const int Flag)
{
    // N = L -1;
    // Definition from Discrete time signal processing. Oppenheim et al., 1999
    const double PI = 3.14159265358979323846;

    switch (Flag)
    {
        case 0 :
                // No window. Apply 1
                // N = L -1; so we just do < instead of <=.
                for (int i = 0; i < N; i++)
                {
                    pWindow[i] = 1.;
                }
                break;
        case 1 :
                // Hamming window.
                for (int i = 0; i < N; i++)
                {
                    pWindow[i] = 0.53836 - 0.46164*cos(2*PI*i/(N-1)); // modified by GR, Dec 2022, because non-symmetrical function
                }
                break;
        case 2 :
                // Hann window.
                for (int i = 0; i < N; i++)
                {
                    pWindow[i] = 0.5 - 0.5*cos(2*PI*i/(N-1)); // modified by GR, Dec 2022, because non-symmetrical function
                }
                break;
        case 3 :
                // for transverse oscillations with two Hamming windows.
            int subaperture_size = 1;
            if (N/3>1) subaperture_size = int(N/3);
                for (int i = 0; i < subaperture_size; i++)
                {
                    pWindow[i] = 0.53836 - 0.46164*cos(2*PI*i/(subaperture_size-1));
                }
                for (int i = N-subaperture_size; i < N; i++)
                {
                    pWindow[i] = 0.53836 - 0.46164*cos(2*PI*(i-(N-subaperture_size))/(subaperture_size-1));
                }
                break;
    }
    return NULL;
}
