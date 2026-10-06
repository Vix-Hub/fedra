//-- Author :  Valeri Tioukov 8/05/2008

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Track momentum estimation algorithms                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TVectorF.h"
#include "TVector.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "EdbLog.h"
#include "EdbPhys.h"
#include "EdbAffine.h"
#include "EdbMomentumEstimator.h"
#include "TFile.h"  
#include <fstream>
#include "TFitResult.h"
#include <TH1F.h>
#include <TH2F.h>

#include <tuple>   
#include <vector>

#include "TRandom3.h"



ClassImp(EdbMomentumEstimator);

using namespace TMath;

//________________________________________________________________________________________
EdbMomentumEstimator::EdbMomentumEstimator()
{
  eAlg   = 0; // default algorithm 
  eF1    = 0;
  eF1X   = 0;
  eF1Y   = 0;
  eG     = 0;
  eGX    = 0;
  eGY    = 0;
  eGA    = 0;
  eGAX   = 0;
  eGAY   = 0;
  eVerbose=0;
  eMinEntr = 1;

  eDTxErrorFun=TF1("dTxError","pol4");
  eDTyErrorFun=TF1("dTyError","pol4");
  eDTsErrorFun=TF1("dTsError","pol4");
  SetParPMS_Mag();
}

//________________________________________________________________________________________
EdbMomentumEstimator::~EdbMomentumEstimator()
{
  SafeDelete(eF1);
  SafeDelete(eF1X);
  SafeDelete(eF1Y);
  SafeDelete(eG);
  SafeDelete(eGX);
  SafeDelete(eGY);
}

//________________________________________________________________________________________
void EdbMomentumEstimator::Set0()
{
  eStatus=-1;
  ePx=ePy=-99;
  eDPx=eDPy=-99;
  ePXmin=ePXmax=ePYmin=ePYmax=-99;  
  eP=eDP=ePmin=ePmax = -99;
}

//________________________________________________________________________________________
void EdbMomentumEstimator::SetParPMS_Mag()
{
  // set the default values for parameters used in PMS_Mag
  eX0 = 4676;       //eX0 = 3504         eX0 = 3521 micron per tung simul      eX0 =  4566 micron per tung+emul simul      eX0 = 4676 micron per tung+emul dati       

  eDTsErrorFun.SetParameters(0.0021, 0.0054,0,0,0);                  
  eDTxErrorFun.SetParameters(0.0021, 0.0093,0,0,0);
  eDTyErrorFun.SetParameters(0.0021, 0.0   ,0,0,0);
}
//________________________________________________________________________________________
void EdbMomentumEstimator::SetParPMS_Mag(Int_t type, Int_t parNumber, Double_t parvalue)
{
  if (type==0) eDTxErrorFun.SetParameter(parNumber,parvalue);
  if (type==1) eDTyErrorFun.SetParameter(parNumber,parvalue);
  if (type==2) eDTsErrorFun.SetParameter(parNumber,parvalue);
}
//________________________________________________________________________________________
void EdbMomentumEstimator::Print()                                
{
  printf("\nEdbMomentumEstimator:\n");
  printf("Algorithm: %s\n", AlgStr(eAlg).Data());
  printf("eX0 = %f \n", eX0);
  printf("eDTxErrorFun parameters:");eDTxErrorFun.Print();
  printf("eDTyErrorFun parameters:");eDTyErrorFun.Print();
  printf("eDTsErrorFun parameters:");eDTsErrorFun.Print();
}

//________________________________________________________________________________________
TString EdbMomentumEstimator::AlgStr(int alg)
{
  char *str = new char[256];
  switch(alg) {
    case 0:
      sprintf(str,"%d (PMSang)        Version rewised by VT 13/05/2008 (see PMSang_base) and further modified by Magali at the end of 2008", alg); break;
    case 1:
      sprintf(str,"%d (PMSang_base)   Version rewised by VT 13/05/2008",alg); break;
    case 2:
      sprintf(str,"%d (PMSang_base_A) Version revised by Andrea Russo 13/03/2009 based on PMSang_base() by VT",alg); break;
    case 3:
      sprintf(str,"%d (PMScoordinate) Momentum estimation by coordinate method A.Russo 2010",alg); break;
    case 4:
      sprintf(str,"%d (PMSang_corr)  PMS_ang corrected by ASh",alg); break;
  }
  TString s(str);
  return s;
}

//________________________________________________________________________________________
float EdbMomentumEstimator::PMS(EdbTrackP &tr)
{
  // according to eAlg the algorithm will be selected to estimate the momentum
  
  eTrack.Clear();
  eTrack.Copy(tr);

  Set0();

  switch(eAlg){
    case 0: return PMSang(eTrack);
    
    case 1: 
      eStatus = PMSang_base(eTrack);
      if(eStatus>0) return eP;
      else return -100.;     // todo

    case 2: 
      eStatus = PMSang_base_A(eTrack);
      if(eStatus>0) return eP;
      else return -100.;     // todo

  //case 3: return PMScoordinate(eTrack);   //ho commentato questo quando ho modificato PMSccordinate in modo che restituisse l'arrey result
  /*case 3: {
    std::vector<float> result = PMScoordinate(eTrack);
    return result[0];  // oppure result[1] se ti serve l'altro valore
}*/
  case 4: return PMSang_corr(eTrack);
  }
  return -100;
}

double CalcTheta(const EdbSegP &s1, const EdbSegP &s2){
  double S1 =s1.TX()*s1.TX()+s1.TY()*s1.TY()+1;
  double S2 =s2.TX()*s2.TX()+s2.TY()*s2.TY()+1;
  double S12=s1.TX()*s2.TX()+s1.TY()*s2.TY()+1;
  double cosTheta2=S12*S12/(S1*S2);
  return 1-cosTheta2;
}

//________________________________________________________________________________________
float EdbMomentumEstimator::PMSang_corr(EdbTrackP &tr)
{
  ////  !Corrections by ASh!
  // Version rewised by VT 13/05/2008 (see PMSang_base) and further modified by Magali at the end of 2008
  // Momentum estimation by multiple scattering (Annecy implementation Oct-2007)
  //
  // Input: tr  - can be modified by the function
  //
  // calculate momentum in transverse and in longitudinal projections using the different 
  // measurements errors parametrisation
   
  int nseg = tr.N();
  int npl=tr.Npl();
  if(nseg<2)   { Log(1,"PMSang","Warning! nseg<2 (%d)- impossible estimate momentum!",nseg);             return -99;}
  if(npl<nseg) { Log(1,"PMSang","Warning! npl<nseg (%d, %d) - use track.SetCounters() first",npl,nseg);  return -99;}
  int plmax = Max( tr.GetSegmentFirst()->PID(), tr.GetSegmentLast()->PID() ) + 1;
  if(plmax<1||plmax>1000)   { Log(1,"PMSang","Warning! plmax = %d - correct the segments PID's!",plmax); return -99;}

  float xmean0,ymean0,zmean0,txmean0,tymean0,wmean0;
  FitTrackLine(tr,xmean0,ymean0,zmean0,txmean0,tymean0,wmean0);    // calculate mean track parameters
  float tmean=Sqrt(txmean0*txmean0+ tymean0*tymean0);

  // -- start calcul --
  const int size = npl;       // vectors size

  TVectorF da(size);
  TArrayI  nentr(size);

  EdbSegP *s1,*s2;
  for(int i1=0; i1<nseg-1; i1++){
        s1 = tr.GetSegment(i1);
        if(!s1) continue;
        
        for(int i2=i1+1; i2<nseg; i2++){
            s2 = tr.GetSegment(i2);
            if(!s2) continue;
            int icell = Abs(s2->PID()-s1->PID());
            da[icell-1] += CalcTheta(*s1,*s2);
            nentr[icell-1] +=1;
        }
    }
 
  float Zcorr = Sqrt(1+tmean*tmean);  // correction due to non-zero track angle and crossed lead thickness

  int max3D=0;                                  // maximum value for the function fit
  TVectorF vind3d(size), errvind3d(size);
  TVectorF errda(size);
  int ist=0;                          // use the counter for case of missing cells 
  for(int i=0; i<size; i++) 
    {
      if( nentr[i] >= eMinEntr/2 && Abs(da[i])<0.1 ) 
        {
          vind3d[ist]    = i+1;                            // x-coord is defined as the number of cells
          errvind3d[ist] = .25;
          da[ist]    = Sqrt( da[i]/(2*nentr[i]) );
          errda[ist] = da[ist]/Sqrt(4*nentr[i]);//CellWeight(npl,i+1)); 
          ist++;
          max3D=i+1;
        }
    }

  float dt = GetDTs(tmean);
  dt*=dt;
  
  float x0    = eX0/(1000*Zcorr);
  float chi2_3D =0;
  
  SafeDelete(eF1);
  SafeDelete(eG);

  eF1 = MCSErrorFunction("eF1",x0,dt);     eF1->SetRange(0,Min(14,max3D));
  eF1->SetParameter(0,2000.);                             // starting value for momentum in GeV

  if (max3D>0)
    {
      eG=new TGraphErrors(vind3d,da,errvind3d,errda);
      if(eG->Fit("eF1","QR")!=0)return -99;
      eP= 0.001*Abs(eF1->GetParameter(0));
      eDP=0.001*eF1->GetParError(0);
      EstimateMomentumError( eP, npl, tmean, ePmin, ePmax );
      chi2_3D = eF1->GetChisquare()/eF1->GetNDF();
      if (eVerbose) printf("P3D=%7.2f GeV ; 90%%C.L. range = [%6.2f : %6.2f] ; chi2_3D %6.2f\n", eP, ePmin, ePmax,chi2_3D);
    }
  return eP;
    // if(eP>50 || eP==2)eP=-99;

  // float ptrue=eP;

  //-----------------------------TO TEST with MC studies------------------------------------
  //  float wx = 1./eDPx/eDPx;
  //  float wy = 1./eDPy/eDPy;
  //  float ptest  = (ePx*wx + ePy*wy)/(wx+wy);  
  //----------------------------------------------------------------------------------------

  // return ptrue;
}
//________________________________________________________________________________________
float EdbMomentumEstimator::PMSang(EdbTrackP &tr)
{
  // Version rewised by VT 13/05/2008 (see PMSang_base) and further modified by Magali at the end of 2008
  //
  // Momentum estimation by multiple scattering (Annecy implementation Oct-2007)
  //
  // Input: tr  - can be modified by the function
  //
  // calculate momentum in transverse and in longitudinal projections using the different 
  // measurements errors parametrisation
   
  int nseg = tr.N();
  int npl=tr.Npl();
  if(nseg<2)   { Log(1,"PMSang","Warning! nseg<2 (%d)- impossible estimate momentum!",nseg);             return -99;}
  if(npl<nseg) { Log(1,"PMSang","Warning! npl<nseg (%d, %d) - use track.SetCounters() first",npl,nseg);  return -99;}
  int plmax = Max( tr.GetSegmentFirst()->PID(), tr.GetSegmentLast()->PID() ) + 1;
  if(plmax<1||plmax>1000)   { Log(1,"PMSang","Warning! plmax = %d - correct the segments PID's!",plmax); return -99;}

 
  float xmean,ymean,zmean,txmean,tymean,wmean;
  float xmean0,ymean0,zmean0,txmean0,tymean0,wmean0;
  FitTrackLine(tr,xmean0,ymean0,zmean0,txmean0,tymean0,wmean0);    // calculate mean track parameters
  float tmean=Sqrt(txmean0*txmean0+ tymean0*tymean0);
  
  EdbSegP *aas;
  float sigmax=0, sigmay=0;
  for (int i =0;i<tr.N();i++)
    {
      aas=tr.GetSegment(i); 
      sigmax+=(txmean0-aas->TX())*(txmean0-aas->TX());   
      sigmay+=(tymean0-aas->TY())*(tymean0-aas->TY());   
    } 
  sigmax= Sqrt(sigmax/tr.N());
  sigmay= Sqrt(sigmay/tr.N());
  for (int i =0;i<tr.N();i++)
    {
      aas=tr.GetSegment(i);
      if(Abs(aas->TX()-txmean0)>3*sigmax||Abs(aas->TY()-tymean0)>3*sigmay) { 
      aas->Set(aas->ID(),aas->X(),aas->Y(),0.,0.,aas->W(),aas->Flag());}
    } 

  FitTrackLine(tr,xmean0,ymean0,zmean0,txmean0,tymean0,wmean0);

  //  EdbAffine2D aff;
  //  aff.ShiftX(-xmean0);
  //  aff.ShiftY(-ymean0);
  //  aff.Rotate( -ATan2(txmean0,tymean0) );   // rotate track to get longitudinal as tx, transverse as ty angle
  //  tr.Transform(aff);
  float PHI=atan2(txmean0,tymean0);
  for (int i =0;i<tr.N();i++)
    {
      aas=tr.GetSegment(i);
      float slx=aas->TY()*cos(-PHI)-aas->TX()*sin(-PHI);
      float sly=aas->TX()*cos(-PHI)+ aas->TY()*sin(-PHI);
      aas->Set(aas->ID(),aas->X(),aas->Y(),slx,sly,aas->W(),aas->Flag());
    }
  
  FitTrackLine(tr,xmean,ymean,zmean,txmean,tymean,wmean);    // calculate mean track parameters
  



  // -- start calcul --

  int minentr  = eMinEntr;               // min number of entries in the cell to accept the cell for fitting
  int stepmax  = npl-1; //npl/minentr;     // max step
  const int size     = stepmax+1;       // vectors size

  TVectorF da(size), dax(size), day(size);
  TArrayI  nentr(size), nentrx(size), nentry(size);

  EdbSegP *s1,*s2;
  for(int ist=1; ist<=stepmax; ist++)         // cycle by the step size
    {
      for(int i1=0; i1<nseg-1; i1++)          // cycle by the first seg
	{
	  s1 = tr.GetSegment(i1);
	  if(!s1) continue;
	  for(int i2=i1+1; i2<nseg; i2++)      // cycle by the second seg
	    {
	      s2 = tr.GetSegment(i2);
	      if(!s2) continue;
	      int icell = Abs(s2->PID()-s1->PID());
	      if( icell == ist ) 
		{
		  if (s2->TX()!=0&&s1->TX()!=0)
		    { 
		      dax[icell-1]   += ( (ATan(s2->TX())- ATan(s1->TX())) * (ATan(s2->TX())- ATan(s1->TX())) );
		      nentrx[icell-1]+=1;
		    }
		  if (s2->TY()!=0&&s1->TY()!=0)
		    { 
		      day[icell-1]   += ( (ATan(s2->TY())- ATan(s1->TY())) * (ATan(s2->TY())- ATan(s1->TY())) );
		      nentry[icell-1]+=1;		   
		    }
		  if (s2->TX()!=0&&s1->TX()!=0&&s2->TY()!=0&&s1->TY()!=0)
		    {
		      da[icell-1]   += (( (ATan(s2->TX())- ATan(s1->TX())) * (ATan(s2->TX())- ATan(s1->TX())) ) 
					+ ( (ATan(s2->TY())-ATan(s1->TY())) * (ATan(s2->TY())- ATan(s1->TY())) ));
		      nentr[icell-1] +=1;
		    }
		}
	    }
	}
    }
 
  float Zcorr = Sqrt(1+txmean0*txmean0+tymean0*tymean0);  // correction due to non-zero track angle and crossed lead thickness

  int maxX =0, maxY=0, max3D=0;                                  // maximum value for the function fit
  TVectorF vindx(size), errvindx(size),vindy(size), errvindy(size),vind3d(size), errvind3d(size);
  TVectorF errda(size), errdax(size), errday(size);
  int ist=0,  ist1=0, ist2=0;                          // use the counter for case of missing cells 
  for(int i=0; i<size; i++) 
    {
      if( nentrx[i] >= minentr && Abs(dax[i])<0.1) 
	{
	  vindx[ist]    = i+1;                            // x-coord is defined as the number of cells
	  errvindx[ist] = .25;
	  dax[ist]    = Sqrt( dax[i]/(nentrx[i]*Zcorr) );
	  errdax[ist] = dax[ist]/Sqrt(2*nentrx[i]);//CellWeight(npl,i+1);    //   Sqrt(npl/vind[i]);
	  ist++;
	  maxX=i+1;
	}
      if( nentry[i] >= minentr && Abs(day[i])<0.1) 
	{
	  vindy[ist1]    = i+1;                            // x-coord is defined as the number of cells
	  errvindy[ist1] = .25;
	  day[ist1]    = Sqrt( day[i]/(nentry[i]*Zcorr) );
	  errday[ist1] = day[ist1]/Sqrt(2*nentry[i]);//CellWeight(npl,i+1);
	  ist1++;
	  maxY=i+1;
	}      
      if( nentr[i] >= minentr/2 && Abs(da[i])<0.1 ) 
	{
	  vind3d[ist2]    = i+1;                            // x-coord is defined as the number of cells
	  errvind3d[ist2] = .25;
	  da[ist2]    = Sqrt( da[i]/(2*nentr[i]*Zcorr) );
	  errda[ist2] = da[ist2]/Sqrt(4*nentr[i]);//CellWeight(npl,i+1));	
	  ist2++;
	  max3D=i+1;
	}
    }

  float dt = GetDTs(tmean);// measurements errors parametrization
  dt*=dt;
  float dtx = GetDTx(txmean);// measurements errors parametrization
  dtx*=dtx;
  float dty = GetDTy(tymean);// measurements errors parametrization
  dty*=dty;
  
  float x0    = eX0/1000;      
  float chi2_3D =0;
  float chi2_T =0;
  float chi2_L =0;
  
  SafeDelete(eF1);
  SafeDelete(eF1X);
  SafeDelete(eF1Y);
  SafeDelete(eG);
  SafeDelete(eGX);
  SafeDelete(eGY);

  eF1X = MCSErrorFunction("eF1X",x0,dtx);    eF1X->SetRange(0,Min(14,maxX));
  eF1X->SetParameter(0,2000.);                             // starting value for momentum in GeV
  eF1Y = MCSErrorFunction("eF1Y",x0,dty);    eF1Y->SetRange(0,Min(14,maxY));
  eF1Y->SetParameter(0,2000.);                             // starting value for momentum in GeV
  eF1 = MCSErrorFunction("eF1",x0,dt);     eF1->SetRange(0,Min(14,max3D));
  eF1->SetParameter(0,2000.);                             // starting value for momentum in GeV

  if (max3D>0)
    {
      eG=new TGraphErrors(vind3d,da,errvind3d,errda);
      eG->Fit("eF1","QR");
      eP=1./1000.*Abs(eF1->GetParameter(0));
      eDP=1./1000.*eF1->GetParError(0);
      if (eP>20||eP<0||eP==2) eP=-99;
      EstimateMomentumError( eP, npl, tymean, ePmin, ePmax );
      chi2_3D = eF1->GetChisquare()/eF1->GetNDF();
      if (eVerbose) printf("P3D=%7.2f GeV ; 90%%C.L. range = [%6.2f : %6.2f] ; chi2_3D %6.2f\n", eP, ePmin, ePmax,chi2_3D);
    }
  if (maxX>0)
    {
      eGX=new TGraphErrors(vindx,dax,errvindx,errdax);
      eGX->Fit("eF1X","QR");
      ePx=1./1000.*Abs(eF1X->GetParameter(0));
      eDPx=1./1000.*eF1X->GetParError(0);
      if (ePx>20||ePx<0||ePx==2) ePx=-99;
      EstimateMomentumError( ePx, npl, txmean, ePXmin, ePXmax );
      chi2_L = eF1X->GetChisquare()/eF1X->GetNDF();
      if (eVerbose) printf("PL=%7.2f GeV ; 90%%C.L. range = [%6.2f : %6.2f] ; chi2_L %6.2f \n",ePx,ePXmin, ePXmax,chi2_L);
    }
  if (maxY>0)
    {
      eGY=new TGraphErrors(vindy,day,errvindy,errday);
      eGY->Fit("eF1Y","QR");
      ePy=1./1000.*Abs(eF1Y->GetParameter(0));
      eDPy=1./1000.*eF1Y->GetParError(0);
      if (ePy>20||ePy<0||ePy==2) ePy=-99;
      EstimateMomentumError( ePy, npl, tmean, ePYmin, ePYmax ); 
      chi2_T = eF1Y->GetChisquare()/eF1Y->GetNDF();
      if (eVerbose) printf("PT=%7.2f GeV ; 90%%C.L. range = [%6.2f : %6.2f] ; chi2_T %6.2f\n", ePy, ePYmin, ePYmax,chi2_T);
    }

  float ptrue=eP;
  if (tmean>0.200&&chi2_T<chi2_3D) 
  {
     ptrue = ePy;
     if (eVerbose) printf(" For this track the evolution of the Transverse projection gives the most accurate estimate of the momentum %7.2f GeV ; 90%%C.L. range = [%6.2f : %6.2f] ; chi2_T /DoF %6.2f\n", ePy, ePYmin, ePYmax,chi2_T);
}

  //-----------------------------TO TEST with MC studies------------------------------------
  //  float wx = 1./eDPx/eDPx;
  //  float wy = 1./eDPy/eDPy;
  //  float ptest  = (ePx*wx + ePy*wy)/(wx+wy);  
  //----------------------------------------------------------------------------------------

  return ptrue;
}

//___________________________________________________________________________________________________


/*float EdbMomentumEstimator::PMScoordinate(EdbTrackP &tr)
{
  // Momentum estimation by coordinate method            
  //
  // April 2010
  
  int nseg = tr.N();
  int npl  = tr.Npl();
  
 
  //float xmean,ymean,zmean,txmean,tymean,wmean;
  float xmean0,ymean0,zmean0,txmean0,tymean0,wmean0;
  FitTrackLine(tr,xmean0,ymean0,zmean0,txmean0,tymean0,wmean0);    // calculate mean track parameters
  float tmean=Sqrt(txmean0*txmean0+ tymean0*tymean0);

  float ang = 0.;

  
  //for (int i =0;i<tr.N();i++)
    //{
      //aas=tr.GetSegment(i);
      //float slx=aas->TY()*cos(-PHI)-aas->TX()*sin(-PHI);
      //float sly=aas->TX()*cos(-PHI)+ aas->TY()*sin(-PHI);
      //aas->Set(aas->ID(),aas->X(),aas->Y(),slx,sly,aas->W(),aas->Flag());
    //}
  
  //FitTrackLine(tr,xmean,ymean,zmean,txmean,tymean,wmean);    // calculate mean track parameters
  

  

  //int minentr  = eMinEntr;               // min number of entries in the cell, should not be set smaller than 5
  int nr1,nr2;

  float dx1,dx2,dy1,dy2,DX1,DY1,DX2,DY2,appx,appy;

  TVectorF da(npl), dax(npl), day(npl);
  TArrayI  nentr(npl);

  for(int i=0;i<npl;i++)
    {
      da[i]    = 0;
      dax[i]   = 0;
      day[i]   = 0;
      nentr[i] = 0;
    }

  EdbSegP *s1=0,*s2=0,*s3=0;

  
  for(int i =0;i<=nseg-3;i++)                   // cycle by the first  seg
    {
      s1 = tr.GetSegment(i);
      for(int j =i+1;j<=nseg-2;j++)             // cycle by the second seg
	{
	  s2 = tr.GetSegment(j);
	  for(int k =j+1;k<=nseg-1;k++)         // cycle by the third  seg
	    {
	      s3 = tr.GetSegment(k);
	      
	      nr1 = TMath::Abs(s1->PID()-s2->PID());
	      nr2 = TMath::Abs(s2->PID()-s3->PID());
	      if(nr1!=nr2) continue;            // continue if (s1,s2) and (s2,s3) correspond to different cells

	      dx1 = s2->X()-s1->X();
	      dx2 = s3->X()-s2->X();
	      dy1 = s2->Y()-s1->Y();
	      dy2 = s3->Y()-s2->Y();
	      
	      ang = 0; // ang is the rotation angle to get trasverse and longitudinal coordinates; not yet implemented 
	      DX1 = cos(ang)*dx1+sin(ang)*dy1;
	      DY1 = cos(ang)*dy1-sin(ang)*dx1;
	      
	      DX2 = cos(ang)*dx2+sin(ang)*dy2;
	      DY2 = cos(ang)*dy2-sin(ang)*dx2;

	      appx = (  DX1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DX2  ) * (  DX1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DX2  );
	      appy = (  DY1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DY2  ) * (  DY1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DY2  );


	      dax[nr1] += appx;
	      day[nr1] += appy;
	      //da[nr1] += (  (((s2->X()-s1->X())*(s2->Z()-s3->Z())/(s1->Z()-s2->Z()))-(s3->X()-s2->X()))**2 + (((y[j]-y[i])*(z[j]-z[k])/(z[i]-z[j]))-(y[k]-y[j]))**2 ) /2. ; 
	      da[nr1] += (appx + appy)/2.; 
	      nentr[nr1] += 1;

      }//end cycle 3rd seg

  }//end cycle 2nd seg

    }//end cycle 1st seg
  
  bool IsEmpty=true; 
  for(int i=0;i<npl;i++)
    {
      if(nentr[i]>0)
  {
    IsEmpty=false;
    dax[i] = sqrt(dax[i]/nentr[i]);
    day[i] = sqrt(day[i]/nentr[i]);
    da[i]  = sqrt(da[i]/nentr[i]);
  }
    }
  
  SafeDelete(eF1);
  SafeDelete(eF1X);
  SafeDelete(eF1Y);
  SafeDelete(eG);
  SafeDelete(eGX);
  SafeDelete(eGY);

  if(IsEmpty)return -99;
  eG  = new TGraphErrors();
  eGX = new TGraphErrors();
  eGY = new TGraphErrors();


  int cont = 0;
  for(int i=0;i<npl;i++)
    {
      if(nentr[i]>eMinEntr)
	{
	  eGX->SetPoint(cont,i*1300,dax[i]);
	  eGX->SetPointError(cont,0,dax[i]/Sqrt(nentr[i]));
	  eGY->SetPoint(cont,i*1300,day[i]);
	  eGY->SetPointError(cont,0,day[i]/Sqrt(nentr[i]));
	  eG->SetPoint(cont,i*1300,da[i]);
	  eG->SetPointError(cont,0,da[i]/Sqrt(nentr[i]));
	  cont++;
	}
    }  
  if(cont==0)return -99;
  eF1X = MCSCoordErrorFunction("eF1X",tmean,eX0);
  eF1X->SetParLimits(0,0.0001,100);
  eF1X->SetParLimits(1,0.0001,100);
  eF1X->SetParameter(0,5);                             // starting value for momentum in GeV
  eF1X->SetParameter(1,10);                              // starting value for coordinate error
  
  eF1Y = MCSCoordErrorFunction("eF1Y",tmean,eX0); 
  //eF1Y->SetRange(0,Min(57,maxY));
  eF1Y->SetParLimits(0,0.0001,100);
  eF1Y->SetParLimits(1,0.0001,100);
  eF1Y->SetParameter(0,5);                             // starting value for momentum in GeV
  eF1Y->SetParameter(1,10);                              // starting value for coordinate error
  
  eF1 = MCSCoordErrorFunction("eF1",tmean,eX0);
  //eF1->SetRange(0,Min(57,max3D));
  eF1->SetParLimits(0,0.0,100);
  eF1->SetParLimits(1,0.0,100);
  eF1->SetParameter(0,5);                             // starting value for momentum in GeV
  eF1->SetParameter(1,10);                              // starting value for coordinate error  
   
  const char *fitopt = "MQ"; //MQR
  eG ->Fit(eF1, fitopt);
  eGX->Fit(eF1X,fitopt);
  eGY->Fit(eF1Y,fitopt);


  eP  = 1./sqrt(eF1->GetParameter(0));
  ePx = 1./sqrt(eF1X->GetParameter(0));
  ePy = 1./sqrt(eF1Y->GetParameter(0));
  return eP;
}*/


float EdbMomentumEstimator::PMScoordinate(EdbTrackP &tr, float sigma_res, float X0, int thickness, bool AddSmear)                               
{ 
  gStyle->SetOptFit(1111);                

  TRandom3 rng(12345); // inizializza generatore Mersenne Twister con seed fisso
  float offset_sigma = 0.15; // sigma dello smearing in micron   


  int trackEvt = tr.GetSegmentFirst()->MCEvt();
     
  int nseg = tr.N();
  int npl  = tr.Npl();  
  
 
  float xmean0,ymean0,zmean0,txmean0,tymean0,wmean0;
  FitTrackLine(tr,xmean0,ymean0,zmean0,txmean0,tymean0,wmean0);    // calculate mean track parameters
  float tmean=TMath::Sqrt(txmean0*txmean0+ tymean0*tymean0);
 
  float ang = 0.;
 
 
  //int minentr  = eMinEntr;               // min number of entries in the cell, should not be set smaller than 5
  int nr1,nr2;
 
  float dx1,dx2,dy1,dy2,DX1,DY1,DX2,DY2,appx,appy;
 
  TVectorF da(npl), dax(npl), day(npl);
  TArrayI  nentr(npl);                  //arrey che tiene traccia di tutte le triplette)
  TArrayI nentr_indep(npl);             //arrey che tiene traccia delle sole triplette indipendenti (da usare per errore su eG)

  std::vector<std::pair<int, int>> used_ranges[npl];    //array (lungo npl) in cui ogni el è un vettore di intervalli di scattering (centerPID, endPID). Serve a tracciare le triplette indip già aggiunte per ciascuna cell length
  
 
  for(int i=0;i<npl;i++)      
    {
      da[i]    = 0;
      dax[i]   = 0;
      day[i]   = 0;
      nentr[i] = 0;
      nentr_indep[i] = 0;
    }
 
  EdbSegP *s1=0,*s2=0,*s3=0;
 
  
  for(int i =0;i<=nseg-3;i++)                   // cycle by the first  seg
    {
      s1 = tr.GetSegment(i);
      for(int j =i+1;j<=nseg-2;j++)             // cycle by the second seg
    {
      s2 = tr.GetSegment(j);
      for(int k =j+1;k<=nseg-1;k++)         // cycle by the third  seg
        {
          s3 = tr.GetSegment(k);
          
          nr1 = TMath::Abs(s1->PID()-s2->PID());
          nr2 = TMath::Abs(s2->PID()-s3->PID());
          if(nr1!=nr2) continue;            // continue if (s1,s2) and (s2,s3) correspond to different cells



          // Applica offset gaussiano a X e Y per simulare disallineamento
          if (AddSmear){
          float x1 = s1->X() + rng.Gaus(0, offset_sigma);
          float x2 = s2->X() + rng.Gaus(0, offset_sigma);
          float x3 = s3->X() + rng.Gaus(0, offset_sigma);

          float y1 = s1->Y() + rng.Gaus(0, offset_sigma);
          float y2 = s2->Y() + rng.Gaus(0, offset_sigma);
          float y3 = s3->Y() + rng.Gaus(0, offset_sigma); 
          
          dx1 = x2 - x1;
          dx2 = x3 - x2;
          dy1 = y2 - y1;
          dy2 = y3 - y2; 
          }
          else{                                                              
          dx1 = s2->X()-s1->X();
          dx2 = s3->X()-s2->X();
          dy1 = s2->Y()-s1->Y();
          dy2 = s3->Y()-s2->Y();  
          }  
             
          
          ang = 0; // ang is the rotation angle to get trasverse and longitudinal coordinates; not yet implemented 
          DX1 = cos(ang)*dx1+sin(ang)*dy1;
          DY1 = cos(ang)*dy1-sin(ang)*dx1;
          
          DX2 = cos(ang)*dx2+sin(ang)*dy2;
          DY2 = cos(ang)*dy2-sin(ang)*dx2;
 
          appx = (  DX1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DX2  ) * (  DX1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DX2  );  
          appy = (  DY1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DY2  ) * (  DY1 * ((s2->Z()-s3->Z())/(s1->Z()-s2->Z())) - DY2  );
 
 
          dax[nr1] += appx;
          day[nr1] += appy;
          //da[nr1] += (  (((s2->X()-s1->X())*(s2->Z()-s3->Z())/(s1->Z()-s2->Z()))-(s3->X()-s2->X()))**2 + (((y[j]-y[i])*(z[j]-z[k])/(z[i]-z[j]))-(y[k]-y[j]))**2 ) /2. ; 
          da[nr1] += (appx + appy)/2.; 
          nentr[nr1] += 1;


          /*int pid2 = s2->PID();   
          int pid3 = s3->PID();
          int start_current = std::min(pid2, pid3);    //individio l'inizio dell'intervallo di scattering per la tripletta corrente (seg centrale)
          int end_current = std::max(pid2, pid3);      //individio la fine dell'intervallo di scattering per la tripletta corrente (seg finale)


          //Verifica se la tripletta corrente è indip (Scopo: Contare solo le triplette indip per ciascuna cell length (cioè per ciascun valore di nr1))
          bool isIndependent = true;                    //variabile booleana che presuppone che la tripletta corrente (s1, s2, s3) sia indip
          for (auto& range : used_ranges[nr1]) {   //itera su tutti gli intervalli di scattering (pid2, pid3) delle triplette già contate come indip per questo valore di cell length (nr1); used_ranges[nr1] è un vettore che contiene per ogni tripletta indip la coppia (pid2, pid3)
            int start = range.first;                //estrae start=pid2 e end=pid3 di ciascuna tripletta indip già salvata per confrontarli con quelli della tripletta corrente
            int end = range.second;
            if (!(end_current <= start || start_current >= end)) {     //pid3 <= start || pid2 >= end è la condizione di NON sovrapposizione degli intervalli di scattering
              isIndependent = false;                   //se c'è sovrapposizone il bool isIndependent viene messo a false
              break;
            }
          }

          if (isIndependent) {                  //se la tripletta corrente è indip allora viene registrata in used_ranges[nr1]    
            used_ranges[nr1].emplace_back(start_current, end_current);     //emplace_back aggiunge la coppia (pid2, pid3) della tripletta corrente ala lista di intervalli di scattering per la cel length nr1
            nentr_indep[nr1] += 1;             //incrementa il contatore delle triplette indip per quel valore di cell length
          }*/
 
      }//end cycle 3rd seg (seg k)
 
  }//end cycle 2nd seg (seg j)
 
    }//end cycle 1st seg (seg i)
  
  bool IsEmpty=true; 
  for(int i=0;i<npl;i++)
    {
      if(nentr[i]>0)
  {
    IsEmpty=false;
    dax[i] = TMath::Sqrt(dax[i]/nentr[i]);
    day[i] = TMath::Sqrt(day[i]/nentr[i]);
    da[i]  = TMath::Sqrt(da[i]/nentr[i]);
  }
    }
  
  //SafeDelete(eF1);
  //SafeDelete(eF1X);   
  //SafeDelete(eF1Y);     
  SafeDelete(eG);
  //SafeDelete(eGX);                              
  //SafeDelete(eGY);               
 
if(IsEmpty)return -99;          


  // Verifica il numero di valori non nulli in nentr
  /*int nonZeroNentr = 0;
  for (int i = 0; i < npl; i++) {
    printf("nentr[%d] = %d (ID: %d)\n", i, nentr[i], tr.ID());
   if (nentr[i] > eMinEntr) nonZeroNentr++;
 }
 printf("per ID %d: Numero di elementi > 0 di nentr= %d\n", tr.ID(), nonZeroNentr);


// Se ci sono esattamente 2 valori maggiori di 0, imposta eP, ePx, ePy a -98
 if (nonZeroNentr == 2) return -98;  */


//DEBUG PER TRIPLETTE INDIPENDENTI
/*printf("\n===== Verifica triplette per track ID %d =====\n", tr.ID());
for (int i = 0; i < npl; ++i) {
  if (nentr[i] > 0) {
    printf("Cell length %2d -> Total: %3d, Independent: %3d\n", i, nentr[i], nentr_indep[i]);
    
    // Stampa delle triplette indipendenti
    if (!used_ranges[i].empty()) {
      printf("  Independent triplets (PID2, PID3):\n");
      for (auto& range : used_ranges[i]) {
        printf("    (%d, %d)\n", range.first, range.second);
      }
    }
  }
}*/




  // Crea i grafici solo se ci sono almeno due valori non nulli
  //if (nonZeroNentr > 1) {
  eG  = new TGraphErrors();                 
  //eGX = new TGraphErrors();
  //eGY = new TGraphErrors();      
  
 
  //int t = 1315; //per simulazioni (60 emulsioni russe per brick)   
  int t = 1350; //per dati (57 emulsioni di Nagoya per brick)                          
  int cont = 0;

  bool firstPointInserted = false;   //bool per identificare il primo punto inserito in eG          

  for(int i=0;i<npl;i++)    
    {
      if(nentr[i]>eMinEntr)     
    {
      // Stampa i valori prima di aggiungerli ai grafici
      //printf("Track ID %d - i = %d, z = %d, da = %.5f, dax = %.5f, day = %.5f, nentr = %d\n",
      //tr.ID(), i, i*t, da[i], dax[i], day[i], nentr[i]);

      /*eGX->SetPoint(cont,i*t,dax[i]);    
      eGX->SetPointError(cont,0,dax[i]/TMath::Sqrt(nentr[i]));
      eGY->SetPoint(cont,i*t,day[i]);
      eGY->SetPointError(cont,0,day[i]/TMath::Sqrt(nentr[i]));        */                    
      eG->SetPoint(cont,i*thickness,da[i]);
      eG->SetPointError(cont,0,da[i]/TMath::Sqrt(2 * nentr[i]));      //HO AGGIUNTO UN 2 SOTTO RADICE PER IL CALCOLO DELL'ERRORE!!!!!!  E ANCHE IL CONTATORE DELLE SOLE TRIPLETTE INDIP nentr_indep[i] AL POSTO DI nentr
      cont++;  
       
        
    }
    }  

   //Istruzioni per costruire la distribuzione dell'RMS del primo punto di eG
      // Inserisci solo il primo valore di eG nell’istogramma ma solo per le tracce che hanno almeno due pt su eG
   /*if (cont > 1 && !firstPointInserted) {
      TH1F *h_Temp_da0 = new TH1F("h_da0", "Distribution of first RMS value; RMS[1];N_{tracks}", 2000, 0, 2); 
      
      double x0, y0;     
      eG->GetPoint(0, x0, y0);  // recupera il primo punto del grafico
      h_Temp_da0->Fill(y0);

      // <<< Qui stampi la info utile  
      //std::cout << "Track ID " << tr.ID() << ": primo punto inserito in eG = da[0] = " << da[0] << std::endl;   

      TFile *outFile = new TFile("RMS_punto0.root", "UPDATE");   
      if (outFile && outFile->IsOpen()) {
        TH1F *hOld_da0 = (TH1F*)outFile->Get("h_da0");
        if (hOld_da0) {  
          hOld_da0->Add(h_Temp_da0);
        } else {
          hOld_da0 = (TH1F*)h_Temp_da0->Clone("h_da0");
        }

        outFile->cd();
        hOld_da0->Write("h_da0", TObject::kOverwrite);
        delete hOld_da0;
        outFile->Close();
      } else {
        std::cerr << "Errore nell'apertura del file ROOT per il salvataggio." << std::endl;       
      }

      delete h_Temp_da0;         
      firstPointInserted = true;    
    }    */                                           

  /*std::cout << "[INFO] Differenze tra ordinate successive di eG (Δy = y_{i+1} - y_i):\n";
  for (int i = 0; i < eG->GetN() - 1; ++i) {
    double x1, y1, x2, y2;
    eG->GetPoint(i, x1, y1);
    eG->GetPoint(i+1, x2, y2);
    double dy = y2 - y1;
    std::cout << "Δy[" << i << "] = y[" << i+1 << "] - y[" << i << "] = "
              << y2 << " - " << y1 << " = " << dy << std::endl;
  }*/
    
  
                                                           
if(cont==0)return -99;   //cont = numero di punti nel grafico
if(cont==1)return -98;                                                     
  
 
  // Calcola il valore di eP, ePx, ePy solo se i grafici sono stati creati                                              
  /*eF1X = MCSCoordErrorFunction("eF1X",tmean,X0); 
  eF1X->SetParLimits(0,0.0001,100);
  eF1X->SetParLimits(1,0.0001,100);
  eF1X->SetParameter(0,5);                             // starting value for momentum in GeV
  eF1X->SetParameter(1,10);    */                          // starting value for coordinate error
    
  
  /*eF1Y = MCSCoordErrorFunction("eF1Y",tmean,X0);
  //TF1* eF1Y = new TF1("eF1Y",Form("sqrt(([1])**2+(2./3)*((x*sqrt(1+%f**2))**3)*(0.0136**2)*[0]/%f)",tmean,X0));
  //std::cout << "[DEBUG] DOPO DEL FIT F1_Y\n";   
  //eF1Y->SetRange(0,Min(57,maxY));
  eF1Y->SetParLimits(0,0.0001,100);                                                                                                          
  eF1Y->SetParLimits(1,0.0001,100);                
  eF1Y->SetParameter(0,5);                             // starting value for momentum in GeV
  eF1Y->SetParameter(1,10);   */                           // starting value for coordinate error
  
      
  //eF1 = MCSCoordErrorFunction("eF1",tmean,X0);                   
  if (!eF1) eF1 = MCSCoordErrorFunction("eF1");   // created once, same formula for all tracks
  eF1->SetParError(0, 0.);                        // no memory of the previous track  
  eF1->SetParLimits(0,0.0,100);          
  eF1->SetParameter(0,5);                             // starting value for momentum in GeV
  eF1->SetParameter(1, TMath::Sqrt(6) * sigma_res);                              // starting value for coordinate error        eF1->SetParameter(1, TMath::Sqrt(6) * 0.15);     

  eF1->FixParameter(1, TMath::Sqrt(6) * sigma_res);                                 //HO AGGIUNTO SQRT(6) AL PARAMETRO DI RISOLUZIONE IN POSIZIONE!!!!!!!      eF1->FixParameter(1, TMath::Sqrt(6) * 0.15);             
  //eF1->SetParLimits(1, TMath::Sqrt(6)*0.5*sigma_res, TMath::Sqrt(6)*1.5*sigma_res);                                                                                          
    
   eF1->SetParameter(2, tmean);                       
   eF1->FixParameter(2, tmean);
   eF1->SetParameter(3, X0);                     
   eF1->FixParameter(3, X0);                                                               

                           
  gErrorIgnoreLevel = kError; // sopprime i warning di Minuit                  
  const char *fitopt = "MQS"; //MQR                                                                   
        


  TFitResultPtr eGResult;                                                                                            
  eGResult = eG->Fit(eF1, fitopt);                                         
  int status = eGResult;                                                                       

 
  if (eGResult->IsValid()){
    eP  = 1./TMath::Sqrt(eF1->GetParameter(0));
    //printf("\n TrackEvt = %d TrackID = %d, eP originario = %.5f\n", eP, trackEvt, tr.ID());  
  } else {
    //std::cout << "[WARN] Fit NON valido per traccia ID: " << tr.ID() << std::endl;  
    eP = -10;
     /*
    // Salvataggio grafico eG in caso di fit divergente
    TFile *divFile = new TFile("MCSgraphs_DivergentFit.root", "UPDATE"); // apri in modalità update
    if (divFile && divFile->IsOpen()) {
      eG->SetName(Form("eG_track%d_div", tr.ID()));  
      eG->SetTitle(Form("eG_Track%d Divergent Fit", tr.ID()));
      eG->Write();
      divFile->Close(); */
    /*std::cout << "Grafico eG salvato in MCSgraphs_divergent.root per fit divergente." << std::endl;*/
  //} else {
    /*std::cerr << "Errore nell'apertura del file ROOT per il salvataggio del fit divergente." << std::endl;*/
  //}
  }

 /*  //COMMENTARE DA QUA
// Calcolo per ogni traccia della variabile SOMMA DEGLI SCARTI NORMALIZZATI PER L'ERRORE e del CHI2
float sumScarti = 0.0;
float chi2 = 0.0;

if (cont > 1 && eP < 1e9) {
  for (int i = 0; i < cont; ++i) {
    double x, y;
    eG->GetPoint(i, x, y);
    double ey = eG->GetErrorY(i);

    // Sicurezza: evitare divisione per zero
    if (ey == 0) continue;

    double valore_fit = eF1->Eval(x);
    double scarto = (y - valore_fit) / ey;

    sumScarti += scarto;
    chi2 += scarto * scarto;
  }


  // --- SCRITTURA SU FILE DI TESTO delle info su chi2 e sumScarti per ogni traccia ---
  int trackEvt = tr.GetSegmentFirst()->MCEvt();
  std::ofstream outFileTxt("Info_Chi_Sum.txt", std::ios::app); // "append" mode
  if (outFileTxt.is_open()) {
    outFileTxt << "EventID: " << trackEvt
             << " nseg: " << nseg
             << "   Chi2: " << chi2
             << "   SumScarti: " << sumScarti << std::endl;
    outFileTxt.close();
  } else {
    std::cerr << "Errore nell'apertura del file Info_Chi_Sum.txt per la scrittura." << std::endl;
  }


  // Istogrammi temporanei
  TH1F *hTempChi2   = new TH1F("hTempChi2", "Chi^{2} Distribution; Chi^{2}; N_{tracks}", 2000, 0, 2000);  
  TH2F *hTemp_sumScarti_chi2 = new TH2F("hTemp_sumScarti_chi2
  TH1F *hTempScarti = new TH1F("hTempScarti", "Distribution of Sum of Normalized Residuals;Sum; N_{tracks}", 2000, -1000, 1000);", "Sum of Norm. Residuals vs Chi^{2}; Chi^{2}; Sum", 2000, 0, 2000, 2000, -1000, 1000);

  hTempScarti->Fill(sumScarti);
  hTempChi2->Fill(chi2);
  hTemp_sumScarti_chi2->Fill(chi2, sumScarti);

  TFile *outFile = new TFile("SumScarti.root", "UPDATE");
  if (outFile && outFile->IsOpen()) {
    // Aggiorna hSumScarti
    TH1F *hOldScarti = (TH1F*)outFile->Get("hSumScarti");
    if (hOldScarti) {
      hOldScarti->Add(hTempScarti);
    } else {
      hOldScarti = (TH1F*)hTempScarti->Clone("hSumScarti");
    }

    // Aggiorna hChi2
    TH1F *hOldChi2 = (TH1F*)outFile->Get("hChi2");
    if (hOldChi2) {
      hOldChi2->Add(hTempChi2);
    } else {
      hOldChi2 = (TH1F*)hTempChi2->Clone("hChi2");
    }

    // Aggiorna h_sumScarti_chi2
    TH2F *hOld_sumScarti_chi2 = (TH2F*)outFile->Get("h_sumScarti_chi2");
    if (hOld_sumScarti_chi2) {
      hOld_sumScarti_chi2 ->Add(hTemp_sumScarti_chi2);
    } else {
      hOld_sumScarti_chi2 = (TH2F*)hTemp_sumScarti_chi2->Clone("h_sumScarti_chi2");
    }

    // Scrivi gli istogrammi nel file
    if (hOldScarti) {
      outFile->cd();
      hOldScarti->Write("hSumScarti", TObject::kOverwrite);
      delete hOldScarti;
    }

    if (hOldChi2) {
      outFile->cd();
      hOldChi2->Write("hChi2", TObject::kOverwrite);
      delete hOldChi2;
    }

    if (hOld_sumScarti_chi2) {
      outFile->cd();
      hOld_sumScarti_chi2->Write("h_sumScarti_chi2", TObject::kOverwrite);
      delete hOld_sumScarti_chi2;
    }

    outFile->Close();
  } else {
    std::cerr << "Errore nell'apertura del file ROOT per il salvataggio." << std::endl;
  }

  delete hTempScarti;
  delete hTempChi2;
  delete hTemp_sumScarti_chi2;




  //REFIT NEL CASO IN CUI IL CHI2 SIA > DEL CHI2_LIMITE PER ALPHA=0.05
  int FreeParameter_eF1 = 1;          //numero di parametri liberi delle funz eF1 = 1 (perchè ho fissato il parametro [1])   
  int ndf = cont - FreeParameter_eF1;   //numero di gradi di liberà del chi2 = n.di pt del TGraph - n. parametri liberi di eF1
  

  if (ndf > 1) {    //ndf>1 equivale a dire che la traccia deve avere cont>2, cioè almeno 3 pt sul Tgraph  
    double chi2_critical = TMath::ChisquareQuantile(0.95, ndf);   //chi2 limite definito da una significatvità alpha=0.05
    if (chi2 > chi2_critical) {
      // Stampa su file il warning
      std::ofstream warnFile("FitWarnings.txt", std::ios::app);  // modalità append
      if (warnFile.is_open()) {
        warnFile << "[WARN] Fit sospetto per traccia Evt " << trackEvt
                << ": chi2 = " << chi2
                << " > chi2_critico = " << chi2_critical << std::endl;
      } else {
        std::cerr << "Errore apertura file FitWarnings.txt" << std::endl;  
      }

      int npoints_refit = (cont > 10) ? 10 : std::max(3, cont / 2);  //n. pt x refit = 10 se la traccia ha più di 10 pt sul Tgraph, è = metà dei pt iniziali (ma almeno pari a 3) se la traccia ha 10 o meno pt sul Tgraph

      TGraphErrors* eG_short = new TGraphErrors();    //creazione di un unovo TGraph con soli 10 punti
      for (int i = 0; i < npoints_refit; ++i) {
        double x, y;
        eG->GetPoint(i, x, y);
        double ex = eG->GetErrorX(i);
        double ey = eG->GetErrorY(i);
        eG_short->SetPoint(i, x, y);
        eG_short->SetPointError(i, ex, ey);
      }

      TFitResultPtr refitResult = eG_short->Fit(eF1, fitopt);
      if (refitResult->IsValid()) {
        eP = 1. / sqrt(eF1->GetParameter(0));          // aggiorna eP con il nuovo fit

        // Ricalcolo chi2 dopo il refit
        double chi2_refit = 0.0;
        for (int i = 0; i < npoints_refit; ++i) {
          double x, y;
          eG_short->GetPoint(i, x, y);
          double ey = eG_short->GetErrorY(i);
          if (ey == 0) continue;
          double yfit = eF1->Eval(x);
          double scarto = (y - yfit) / ey;
          chi2_refit += scarto * scarto;
        }

        // Scrive anche chi2_refit ed eP nel file
        if (warnFile.is_open()) {
          warnFile << "       --> Nuovo chi2 dopo refit con " << npoints_refit
                  << " punti: " << chi2_refit << std::endl;
          warnFile << "       --> [INFO] Fit limitato valido. Nuovo eP = " << eP << "\n";
          warnFile.close();                 
        }

        // Salva grafico refittato nel file separato
        TFile *refitFile = new TFile("MCSgraphs_refit.root", "UPDATE");
        if (refitFile && refitFile->IsOpen()) {
          eG_short->SetName(Form("eG_refit_track%d", trackEvt));  //tr.ID()
          eG_short->SetTitle(Form("Refit - Track Evt %d", trackEvt));  //tr.ID()
          eG_short->Write();
          refitFile->Close();
        } else {
          std::cerr << "Errore apertura file per salvataggio refit.\n";
        }
      } else {
        std::cerr << "[WARN] Refit su 10 punti fallito per traccia Evt: " << trackEvt << "\n"; //tr.ID()
      }

      delete eG_short; // Libera memoria
}

}
} */   //A QUA

  //SCOMMENTARE DA QUA 
//TH1D *hTemp_DecreasingY = new TH1D("hTemp_DecreasingY", "Distribution of #Delta y (y_{i+1} < y_{i});#Delta y [#mum];Counts", 1000, -10, 10);

//const float maxDecrease = 0.1; // tolleranza  0.1 micron  
const float minIncrease = 0.1;        

const double alpha = 0.8;                                                         

TGraphErrors* eG_crescente = new TGraphErrors();      
int nPoints = eG->GetN();
    
double x_prev, y_prev;
eG->GetPoint(0, x_prev, y_prev);     
double ex_prev = eG->GetErrorX(0);            
double ey_prev = eG->GetErrorY(0);                 
      
eG_crescente->SetPoint(0, x_prev, y_prev);
eG_crescente->SetPointError(0, ex_prev, ey_prev);

int pt_count = 1;
for (int i = 1; i < nPoints; ++i) {
    double x, y;
    eG->GetPoint(i, x, y);
    double ex = eG->GetErrorX(i);
    double ey = eG->GetErrorY(i);

    //if (y >= y_prev - maxDecrease) {    
    if (y >= y_prev + minIncrease) {    
    //if (y >= y_prev - alpha * ey){
        eG_crescente->SetPoint(pt_count, x, y);
        eG_crescente->SetPointError(pt_count, ex, ey);           
        y_prev = y;
        ++pt_count;
    } else {
      //double deltaY = y - y_prev;
      //hTemp_DecreasingY->Fill(deltaY);
        // Appena un punto non soddisfa la condizione, interrompe il ciclo di riempimento di eG_crescente
        break;
    }
}

//printf("Pt su eG_crescente (prima di iter.): %d\n", pt_count);  //A QUA         


/*lasciare commentato
if (pt_count > 2) {
    TFitResultPtr refitResult = eG_crescente->Fit(eF1, fitopt);
    if (refitResult->IsValid()) {
        eP = 1. / sqrt(eF1->GetParameter(0));  // aggiorna eP con il nuovo fit

        std::cout << "Differenze (residui) tra punti e curva di fit per traccia " << trackEvt << ":\n";
        const double residuoThreshold = 0.2;
        bool refitNeeded = false;

        int nCrescentPoints = eG_crescente->GetN();
        int idx = nCrescentPoints - 1;

        // Controlla partendo dall'ultimo punto e rimuove i punti finali con residuo > soglia
        while (idx >= 0 && eG_crescente->GetN() > 2) {
            double x, y;
            eG_crescente->GetPoint(idx, x, y);
            double y_fit = eF1->Eval(x);
            double residuo = y - y_fit;

            std::cout << "Punto " << idx << ": x = " << x << ", y = " << y << ", y_fit = " << y_fit << ", residuo = " << residuo << "\n";

            if (std::abs(residuo) > residuoThreshold) {
                std::cout << "[INFO] Punto " << idx << " ha residuo > 0.2 e sarà rimosso.\n";
                eG_crescente->RemovePoint(idx);
                refitNeeded = true;
                idx--;  // Passa al punto precedente
            } else {
                break;  // Interrompe il controllo se il punto è accettabile
            }
        }

        // Se sono stati rimossi punti, esegui di nuovo il fit
        if (refitNeeded && eG_crescente->GetN() > 2) {
            std::cout << "[INFO] Eseguo refit dopo rimozione punti anomali...\n";
            refitResult = eG_crescente->Fit(eF1, fitopt);
            if (refitResult->IsValid()) {
                eP = 1. / sqrt(eF1->GetParameter(0));
                std::cout << "[INFO] Refit riuscito dopo rimozione outlier.\n";

                // --- Calcolo e stampa dei nuovi residui ---
                std::cout << "Residui DOPO refit per traccia " << trackEvt << ":\n";
                int newN = eG_crescente->GetN();
                for (int i = 0; i < newN; ++i) {
                    double x, y;
                    eG_crescente->GetPoint(i, x, y);
                    double y_fit = eF1->Eval(x);
                    double residuo = y - y_fit;
                    std::cout << "Punto " << i << ": x = " << x << ", y = " << y << ", y_fit = " << y_fit << ", residuo = " << residuo << "\n";
                }

                // Salva grafico refittato nel file separato
                TFile *refitFile = new TFile("MCSgraphs_crescenti_refit.root", "UPDATE");
                if (refitFile && refitFile->IsOpen()) {
                    eG_crescente->SetName(Form("eG_crescente_track%d", trackEvt));
                    eG_crescente->SetTitle(Form("Refit - Track Evt %d", trackEvt));
                    eG_crescente->Write();
                    refitFile->Close();
                } else {
                    std::cerr << "Errore apertura file per salvataggio refit.\n";
                }

            } else {
                std::cerr << "[WARN] Refit fallito dopo rimozione outlier per traccia Evt: " << trackEvt << "\n";
            }
        }
    } else {
        std::cerr << "[WARN] Refit TGraph crescente fallito per traccia Evt: " << trackEvt << "\n";
    }
}*/

   //SCOMMENTARE DA QUA
/*auto rimuovi_outlier_finali = [&](TGraphErrors* graph, TF1* fitFunc, const double threshold) -> bool {
    bool puntiRimossi = false;
    int idx = graph->GetN() - 1;

    while (idx >= 0 && graph->GetN() > 2) {  
        double x, y;
        graph->GetPoint(idx, x, y);
        double y_fit = fitFunc->Eval(x);
        double residuo = y - y_fit;

        //std::cout << "Punto " << idx << ": x = " << x << ", y = " << y << ", y_fit = " << y_fit << ", residuo = " << residuo << "\n";

        if (std::abs(residuo) > threshold) {
            //std::cout << "[INFO] Punto " << idx << " ha residuo > " << threshold << " e sarà rimosso.\n";
            graph->RemovePoint(idx);
            puntiRimossi = true;
            idx--;
        } else {
            break;
        }
    }

    return puntiRimossi;
}; */

//Versione in cui rimuovo anche outliers intermedi:
auto rimuovi_outlier_finali = [&](TGraphErrors* graph, TF1* fitFunc, const double threshold) -> bool {
    bool puntiRimossi = false;

    for (int idx = graph->GetN() - 1; idx >= 0; --idx) {     

        if (graph->GetN() <= 3) break;         

        double x, y;             
        graph->GetPoint(idx, x, y);
        double y_fit = fitFunc->Eval(x);
        double ey = graph->GetErrorY(idx);   //riga aggiunta per calcolare residuo in unità di sigma
        //double residuo = (y - y_fit) / ey;
        double residuo = y - y_fit;

        if (std::abs(residuo) > threshold) {
            graph->RemovePoint(idx);
            puntiRimossi = true;
        }
    }

    return puntiRimossi;
};


TGraphErrors* eG_crescente_iniziale = nullptr; 
TGraphErrors* eG_crescente_post1 = nullptr;
TGraphErrors* eG_crescente_post2 = nullptr; 


if (eG_crescente->GetN() > 2) {
    TFitResultPtr refitResult = eG_crescente->Fit(eF1, fitopt);   //se non vengono eseguite le iterazioni 1 e 2 l'impulso restituito è quello ottenuto dal grafico eG_crescente
    if (refitResult->IsValid()) {
        eP = 1. / TMath::Sqrt(eF1->GetParameter(0));

        // Salva stato prima della prima iterazione, cioè il grafico crescente iniziale:
        eG_crescente_iniziale = (TGraphErrors*)eG_crescente->Clone("eG_crescente_iniziale");

        const double residuoThreshold_first = 0.2;   
        const double residuoThreshold_second = 0.2;                                                    

        bool primoRefitRiuscito = false;         

        //std::cout << "[INFO] Inizio prima iterazione rimozione outlier...\n";
        bool rimossi_1 = rimuovi_outlier_finali(eG_crescente, eF1, residuoThreshold_first);

        bool skipSecondIteration = false;

        if (rimossi_1) {
            int remainingPoints = eG_crescente->GetN();
            if (remainingPoints > 2) {
                refitResult = eG_crescente->Fit(eF1, fitopt);
                if (refitResult->IsValid()) {
                    eP = 1. / TMath::Sqrt(eF1->GetParameter(0));
                    //std::cout << "[INFO] Primo refit riuscito.\n";
                    primoRefitRiuscito = true;
                    // Salva stato dopo la prima iterazione:
                    eG_crescente_post1 = (TGraphErrors*)eG_crescente->Clone("eG_crescente_post1");
                } else {
                    //std::cerr << "[WARN] Primo refit fallito.\n";
                    skipSecondIteration = true;
                }
            } else {
                //std::cerr << "[WARN] Troppi pochi punti (" << remainingPoints << ") dopo prima rimozione. Refitting saltato.\n";
                //std::cout << "[INFO] Seconda iterazione saltata per numero insufficiente di punti.\n";   
                skipSecondIteration = true;      
            }
        } else { 
            skipSecondIteration = true;
        }            

        if (!skipSecondIteration) {
            //std::cout << "[INFO] Inizio seconda iterazione rimozione outlier...\n";
            bool rimossi_2 = rimuovi_outlier_finali(eG_crescente, eF1, residuoThreshold_second);

            if (rimossi_2) {
                int remainingPoints = eG_crescente->GetN();
                if (remainingPoints > 2) {
                    refitResult = eG_crescente->Fit(eF1, fitopt);
                    if (refitResult->IsValid()) {
                        eP = 1. / TMath::Sqrt(eF1->GetParameter(0));
                        //std::cout << "[INFO] Secondo refit riuscito.\n";

                        // Salva stato dopo la seconda iterazione:
                        eG_crescente_post2 = (TGraphErrors*)eG_crescente->Clone("eG_crescente_post2");
                    } else {
                        //std::cerr << "[WARN] Secondo refit fallito.\n";
                    }
                } else {
                    //std::cerr << "[WARN] Troppi pochi punti (" << remainingPoints << ") dopo seconda rimozione. Refitting saltato.\n";
                }
            }
        }   //A QUA     

        /*   // Stampa residui finali
        std::cout << "Residui FINALI dopo due iterazioni per traccia " << trackEvt << ":\n";
        int finalN = eG_crescente->GetN();
        for (int i = 0; i < finalN; ++i) {
            double x, y;
            eG_crescente->GetPoint(i, x, y);
            double y_fit = eF1->Eval(x);
            double residuo = y - y_fit;
            std::cout << "Punto " << i << ": x = " << x << ", y = " << y << ", y_fit = " << y_fit << ", residuo = " << residuo << "\n";
        }*/

        /* // Salvataggio file
        TFile *refitFile = new TFile("MCSgraphs_crescenti_refit.root", "UPDATE");
        if (refitFile && refitFile->IsOpen()) {
            eG_crescente->SetName(Form("eG_crescente_track%d", trackEvt));
            eG_crescente->SetTitle(Form("Refit - Track Evt %d", trackEvt));
            eG_crescente->Write();
            refitFile->Close();
        } else {
            std::cerr << "Errore apertura file per salvataggio refit.\n";
        }
    } else {
        std::cerr << "[WARN] Fit iniziale fallito per traccia Evt: " << trackEvt << "\n";
    }*/ 

 //SCOMMENTARE DA QUA
// Salvataggio file condizionato  
/*TFile *refitFile = new TFile("MCSgraphs_crescente_refit.root", "UPDATE");      
if (refitFile && refitFile->IsOpen()) {
    if (eG_crescente_post2) {
        //std::cout << "[INFO] Salvando il grafico dopo la seconda iterazione.\n";
        eG_crescente_post2->SetName(Form("eG_crescente_post2_track%d", trackEvt));
        eG_crescente_post2->SetTitle(Form("Refit - Track Evt %d (TrackID %d) - Second Iteration", trackEvt, tr.ID()));
         eG_crescente_post2->Draw("AP"); 
         gPad->Update();
         gPad->GetListOfPrimitives()->Write("stats", TObject::kOverwrite);
        eG_crescente_post2->Write();
    } else if (eG_crescente_post1) {
        //std::cout << "[INFO] Salvando il grafico dopo la prima iterazione.\n";
        eG_crescente_post1->SetName(Form("eG_crescente_post1_track%d", trackEvt));
        eG_crescente_post1->SetTitle(Form("Refit - Track Evt %d (TrackID %d) - First Iteration", trackEvt, tr.ID()));
         eG_crescente_post1->Draw("AP");
        gPad->Update();
        gPad->GetListOfPrimitives()->Write("stats", TObject::kOverwrite);
        eG_crescente_post1->Write();
    } else {              
        //std::cout << "[INFO] Salvando il grafico crescente iniziale, senza iterazioni.\n";
        eG_crescente_iniziale->SetName(Form("eG_crescente_track%d", trackEvt));
        eG_crescente_iniziale->SetTitle(Form("Increasing Graph (Before Iter.) - Track Evt %d (TrackID %d)", trackEvt, tr.ID()));
        eG_crescente_iniziale->Draw("AP");
        gPad->Update();
        gPad->GetListOfPrimitives()->Write("stats", TObject::kOverwrite);   
        eG_crescente_iniziale->Write();
    }
    refitFile->Close();
} else {       
    std::cerr << "Errore apertura file per salvataggio refit.\n";
}        */                                                                                              
   

    }
      delete eG_crescente_iniziale; // Libera memoria
      delete eG_crescente_post1;
      delete eG_crescente_post2;
      delete eG_crescente;  
  
  
  }   

           
else {

    //------------------------------------------------------------
    // ⬇ NUOVO BLOCCO: ITERAZIONI ANCHE QUANDO eG_crescente NON ESISTE
    // eG_clone = copia completa del grafico originale
    TGraphErrors* eG_clone = (TGraphErrors*)eG->Clone("eG_clone");

    TGraphErrors* eG_clone_iniziale = nullptr;
    TGraphErrors* eG_clone_post1     = nullptr;
    TGraphErrors* eG_clone_post2     = nullptr;         
     
    // Fit iniziale su clone
    TFitResultPtr refitResult = eG_clone->Fit(eF1, fitopt);
    if (refitResult->IsValid()) {

        eP = 1. / TMath::Sqrt(eF1->GetParameter(0));  

        eG_clone_iniziale = (TGraphErrors*)eG_clone->Clone("eG_clone_iniziale");

        const double residuoThreshold_first = 0.2;
        const double residuoThreshold_second = 0.2;                                               
     
        bool rimossi_1 = rimuovi_outlier_finali(eG_clone, eF1, residuoThreshold_first);           
        bool skipSecondIteration = false;

        if (rimossi_1 && eG_clone->GetN() > 2) {
            refitResult = eG_clone->Fit(eF1, fitopt);
            if (refitResult->IsValid()) {
                eP = 1. / TMath::Sqrt(eF1->GetParameter(0));
                eG_clone_post1 = (TGraphErrors*)eG_clone->Clone("eG_clone_post1");
            } else skipSecondIteration = true;
        } else skipSecondIteration = true;

        if (!skipSecondIteration) {
            bool rimossi_2 = rimuovi_outlier_finali(eG_clone, eF1, residuoThreshold_second);
            if (rimossi_2 && eG_clone->GetN() > 2) {
                refitResult = eG_clone->Fit(eF1, fitopt);
                if (refitResult->IsValid()) {
                    eP = 1. / TMath::Sqrt(eF1->GetParameter(0));
                    eG_clone_post2 = (TGraphErrors*)eG_clone->Clone("eG_clone_post2");
                }
            }           
        }

        // SALVATAGGIO FILE
        
        /*TFile *refitFile = new TFile("MCSgraphs_crescente_refit.root", "UPDATE");
        if (refitFile && refitFile->IsOpen()) {

            if (eG_clone_post2) {
                eG_clone_post2->SetName(Form("eG_clone_post2_track%d", trackEvt));
                eG_clone_post2->SetTitle(Form("eG_clone_post2 - Track Evt %d (TrackID %d)", trackEvt, tr.ID()));
                eG_clone_post2->Draw("AP");
		gPad->Update();
		gPad->GetListOfPrimitives()->Write("stats", TObject::kOverwrite);
                eG_clone_post2->Write();
            }
            else if (eG_clone_post1) {
                eG_clone_post1->SetName(Form("eG_clone_post1_track%d", trackEvt));
                eG_clone_post1->SetTitle(Form("eG_clone_post1 - Track Evt %d (TrackID %d)", trackEvt, tr.ID()));
                eG_clone_post1->Draw("AP");
		gPad->Update();
		gPad->GetListOfPrimitives()->Write("stats", TObject::kOverwrite);
                eG_clone_post1->Write();           
            }

            refitFile->Close();
        }    */  

        delete eG_clone_iniziale;
        delete eG_clone_post1;
        delete eG_clone_post2;
    }
}
    

  
  /*TFitResultPtr eGXResultX = eGX->Fit(eF1X,fitopt);
  status = eGXResultX;
  if (status==0) {  */
    /*std::cout << "Valid Fit" << std::endl;*/
    /*ePx  = 1./TMath::Sqrt(eF1->GetParameter(0));
  } else {  */
    /*std::cout << "Divergent Fit" << std::endl;*/
    /*ePx = -10;
  }  */    

  /*TFitResultPtr eGYResultY = eGY->Fit(eF1Y,fitopt);                                   
  status = eGYResultY;
  if (status==0) {   */
    /*std::cout << "Valid Fit" << std::endl;*/
    /*ePy  = 1./TMath::Sqrt(eF1->GetParameter(0));
  } else {  */
    /*std::cout << "Divergent Fit" << std::endl;*/   
    /*ePy = -10;    
  }   */
//}

 
  /*eP  = 1./sqrt(eF1->GetParameter(0));
  ePx = 1./sqrt(eF1X->GetParameter(0));   
  ePy = 1./sqrt(eF1Y->GetParameter(0));*/



 // Creazione e salvataggio dei grafici su un file ROOT
/*TFile *outFile = new TFile("MCSgraphs.root", "UPDATE"); //apre in modalità update
if (outFile && outFile->IsOpen()) {
  eG->SetName(Form("eG_track%d", tr.GetSegmentFirst()->MCEvt() ));   //eG->SetName(Form("eG_track%d", tr.ID()));  
  eG->SetTitle(Form("Track Evt %d (TrackID %d)", tr.GetSegmentFirst()->MCEvt(), tr.ID() ));   //eG->SetTitle(Form("eG_Track%d", tr.ID()));

  //eGX->SetName(Form("eGX_track%d", tr.ID()));
  // eGY->SetName(Form("eGY_track%d", tr.ID())); 
  
  eG->Draw("AP");
  gPad->Update();
  gPad->GetListOfPrimitives()->Write("stats", TObject::kOverwrite);   
      
  eG->Write();
  //eGX->Write();                   
  //eGY->Write();               
             
  outFile->Close();   */
  /*std::cout << "Grafici salvati su MCSgraphs.root" << std::endl;*/
//} else {      
  /*std::cerr << "Errore nell'apertura del file ROOT per il salvataggio." << std::endl;*/    
//}                                            
                            

return eP;                              
}             

//________________________________________________________________________________________
float EdbMomentumEstimator::CellWeight(int npl, int m)
{
  //--------------------------- TO BE IMPLEMENTED-----------------------------------------

  // npl - number of plates, m - the cell thickness in plates
  // return the statistical weight of the cell

  //return  Sqrt(npl/m);  // the simpliest estimation no shift, no correlations

  return 2*Sqrt( npl/m + 1./m/m*( npl*(m-1) - m*(m-1)/2.) );
  // return 1;
}

//________________________________________________________________________________________
TF1 *EdbMomentumEstimator::MCSErrorFunction(const char *name, float x0, float dtx)
{
  //        dtx - the plane angle measurement error
  // return the function of the expected angular deviation vs range
  //
  // use the Highland-Lynch-Dahl formula for theta_rms_plane = 13.6 MeV/bcp*z*sqrt(x/x0)*(1+0.038*log(x/x0))  (PDG)
  // so the expected measured angle is sqrt( theta_rms_plane**2 + dtx**2)
  //
  // The constant term im the scattering formula is not 13.6 but 14.64, which
  // is the right reevaluated number, due to a calculation with the moliere
  // distribution. 13.6 is an approximation. See Geant3 or 4 references for more explanations.???????
  //
  // err(x) = sqrt(k*x*(1+0.038*log(x/x0))/p**2 + dtx)

  //  float k   = 14.64*14.64/x0;
  // 14.64*14.64/1000/1000 = 0.0002143296  - we need p in GeV 
  // 13.6*13.6/1000/1000   = 0.0001849599  - we need p in GeV


  return new TF1(name,Form("sqrt(214.3296*x/%f*((1+0.038*log(x/(%f)))**2)/([0])**2+%e)",x0,x0,dtx));
  // return new TF1(name,Form("sqrt(184.9599*x/%f*((1+0.038*log(x/(%f)))**2)/([0])**2+%e)",x0,x0,dtx));             

  //P is returned in MeV by this function for more convinience, but given in GeV as output.
}

//________________________________________________________________________________________

/*TF1 *EdbMomentumEstimator::MCSCoordErrorFunction(const char *name, float tmean,float x0)
{
  // return the function of the expected position deviation as function of range
  //
  // use the Highland-Lynch-Dahl formula for theta_rms_plane = 13.6 MeV/bcp*z*sqrt(x/x0)  (PDG)
  // so the expected measured position deviation is (1/sqrt(3))*theta_rms_plane + measurement error
  // log term in HLD formula is neglected  
       

  return new TF1(name,Form("sqrt(([1])**2+(2./3)*((x*sqrt(1+%f**2))**3)*(0.0136**2)*[0]/%f)",tmean,x0));                 
  
  //P is returned in GeV
}    */

TF1 *EdbMomentumEstimator::MCSCoordErrorFunction(const char *name)
{
  // return the function of the expected position deviation as function of range
  //
  // use the Highland-Lynch-Dahl formula for theta_rms_plane = 13.6 MeV/bcp*z*sqrt(x/x0)  (PDG)
  // so the expected measured position deviation is (1/sqrt(3))*theta_rms_plane + measurement error
  // log term in HLD formula is neglected  
       

  return new TF1(name,Form("sqrt(([1])**2+(2./3)*((x*sqrt(1+[2]**2))**3)*(0.0136**2)*[0]/[3])"), 0.0, 1000000.0);                 
  
  //P is returned in GeV
} 



/*TF1 *EdbMomentumEstimator::MCSCoordErrorFunction(const char *name, float tmean, float x0)
{
  // Expected position deviation as function of range (Highland–Lynch–Dahl)
  // theta_rms_plane = 13.6 MeV/bcp*z*sqrt(x/x0)
  // position deviation = (1/sqrt(3))*theta_rms_plane + measurement error             
  // log term neglected

  TString expr = Form(
    "sqrt( pow([1],2) + (2./3) * pow( x * pow(1 + %g*%g, 0.5) , 3) * pow(0.0136,2) * [0] / %g )",
    (double)tmean, (double)tmean, (double)x0
  );

  printf("DEBUG: MCSCoordErrorFunction formula = %s\n", expr.Data());

  return new TF1(name, expr.Data());
}*/

/*TF1 *EdbMomentumEstimator::MCSCoordErrorFunction(const char *name, float tmean,float x0)
{
  // return the function of the expected position deviation as function of range
  //
  // use the Highland-Lynch-Dahl formula for theta_rms_plane = 13.6 MeV/bcp*z*sqrt(x/x0)  (PDG)
  // so the expected measured position deviation is (1/sqrt(3))*theta_rms_plane + measurement error
  // log term in HLD formula is neglected  
       

  return new TF1(name,Form("TMath::Sqrt(([1])**2+(2./3)*((x*TMath::Sqrt(1+%f**2))**3)*(0.0136**2)*[0]/%f)",tmean,x0));                 
  
  //P is returned in GeV                                                 
} */       

/*TF1 *EdbMomentumEstimator::MCSCoordErrorFunction(const char *name, float tmean, float x0)
{
    TString expr = Form(
        "sqrt( [1]^2 + (2./3) * ( x * sqrt(1 + (%g)^2) )^3 * (0.0136)^2 * [0] / %g )",
        tmean, x0
    );

    printf("DEBUG: MCSCoordErrorFunction formula = %s\n", expr.Data());

    return new TF1(name, expr.Data());
} */                                                             
                                                                                                          


//______________________________________________________________________________________
void EdbMomentumEstimator::EstimateMomentumError(float P, int npl, float ang, float &pmin, float &pmax)
{
  float pinv=1./P;
  float  DP=Mat(P, npl, ang );
  float pinvmin=pinv*(1-DP*1.64);
  float pinvmax=pinv*(1+DP*1.64);
  pmin=(1./pinvmax);   //90%CL minimum momentum
  pmax=(1./pinvmin);   //90%CL maximum momentum
  if (P>1000.) pmax=10000000.;
}

//______________________________________________________________________________________
double EdbMomentumEstimator::Mat(float P, int npl, float ang)
{
  // These parametrisations at low and large angles are parametrised with MC
  // See Magali's thesis for more informations
  double DP=0.;
// new parameterisition as published
	if(Abs(ang)<0.2)   DP = ((0.397+0.019*P)/Sqrt(npl) + (0.176+0.042*P) + (-0.014-0.003*P)*Sqrt(npl));
	if(Abs(ang)>=0.2)  DP = ((1.400-0.022*P)/Sqrt(npl) + (-0.040+0.051*P) + (0.003-0.004*P)*Sqrt(npl));
  if (DP>0.80) DP=0.80;
  return DP;
}

//________________________________________________________________________________________
void EdbMomentumEstimator::DrawPlots(TCanvas *c1)
{
  // example of the plots available after PMSang
  gStyle->SetOptFit(11111);
  if(c1==NULL) c1 = new TCanvas("tf","track momentum estimation",800,600);
  c1->Divide(3,2);

  if(eAlg<2||eAlg==3)
    {
      if(eGX) {
	c1->cd(1);
	TGraphErrors *gx = new TGraphErrors(*eGX);

	if(eAlg<2)
	  {
	    gx->Draw("ALPR");
	    gx->SetTitle("Theta vs cell (transverse component)");
	    TF1 *fxmin = new TF1(*(eF1X));
	    fxmin->SetLineColor(kBlue);
	    fxmin->SetParameter(0,ePXmin);
	    fxmin->Draw("same");
	    TF1 *fxmax = new TF1(*(eF1X));
	    fxmax->SetLineColor(kBlue);
	    fxmax->SetParameter(0,ePXmax);
	    fxmax->Draw("same");
	  }
	else
	  {
	    gx->Draw("A*");
	    gx->SetTitle("Position vs cell (X component)");
	  }
	
	
      }
      
      if(eGY) {
	c1->cd(2);
	TGraphErrors *gy = new TGraphErrors(*eGY);

	if(eAlg<2)
	  {
	    gy->Draw("ALPR");
	    gy->SetTitle("Theta vs cell (longitudinal component)");
	    TF1 *fymin = new TF1(*(eF1Y));
	    fymin->SetLineColor(kBlue);
	    fymin->SetParameter(0,ePYmin);
	    fymin->Draw("same");
	    TF1 *fymax = new TF1(*(eF1Y));
	    fymax->SetLineColor(kBlue);
	    fymax->SetParameter(0,ePYmax);
	    fymax->Draw("same");
	  }
	else
	  {
	    gy->Draw("A*");
	    gy->SetTitle("Position vs cell (Y component)");
	  }
	
	
      }
      
      if(eG) {
	c1->cd(3);
	TGraphErrors *g = new TGraphErrors(*eG);

	if(eAlg<2)
	  {
	    g->Draw("ALPR");
	    g->SetTitle("Theta vs cell (3D)");
	    //g->Print();
	    TF1 *fmin = new TF1(*(eF1));
	    fmin->SetLineColor(kBlue);
	    fmin->SetParameter(0,ePmin);
	    fmin->Draw("same");
	    TF1 *fmax = new TF1(*(eF1));
	    fmax->SetLineColor(kBlue);
	    fmax->SetParameter(0,ePmax);
	    fmax->Draw("same");
	  }
	else
	  {
	    g->Draw("A*");
	    g->SetTitle("Position vs cell (3D)");
	  }
 
	
      }
    }

  if(eAlg==2)
    {
      if(eGAX) {
	c1->cd(1);
	TGraphAsymmErrors *gx = new TGraphAsymmErrors(*eGAX);
	gx->SetTitle("Theta vs cell (longitudinal component)");
	gx->Draw("ALPR");
	TF1 *fxmin = new TF1(*(eF1X));
	fxmin->SetLineColor(kBlue);
	fxmin->SetParameter(0,ePXmin);
	fxmin->Draw("same");
	TF1 *fxmax = new TF1(*(eF1X));
	fxmax->SetLineColor(kBlue);
	fxmax->SetParameter(0,ePXmax);
	fxmax->Draw("same");
      }
      
      if(eGAY) {
	c1->cd(2);
	TGraphAsymmErrors *gy = new TGraphAsymmErrors(*eGAY);
	gy->SetTitle("Theta vs cell (transverse component)");
	gy->Draw("ALPR");
	TF1 *fymin = new TF1(*(eF1Y));
	fymin->SetLineColor(kBlue);
	fymin->SetParameter(0,ePYmin);
	fymin->Draw("same");
	TF1 *fymax = new TF1(*(eF1Y));
	fymax->SetLineColor(kBlue);
	fymax->SetParameter(0,ePYmax);
	fymax->Draw("same");
      }
      
      if(eGA) {
	c1->cd(3);
	TGraphAsymmErrors *g = new TGraphAsymmErrors(*eGA);
	g->SetTitle("Theta vs cell (3D)");
	g->Draw("ALPR");
	//g->Print();
	TF1 *fmin = new TF1(*(eF1));
	fmin->SetLineColor(kBlue);
	fmin->SetParameter(0,ePmin);
	fmin->Draw("same");
	TF1 *fmax = new TF1(*(eF1));
	fmax->SetLineColor(kBlue);
	fmax->SetParameter(0,ePmax);
	fmax->Draw("same");
      }
    }


  
      int nseg = eTrack.N();
      if(nseg) {
	c1->cd(4);
    TGraph2D *gxyz = new TGraph2D(nseg);
    for(int i=0; i<nseg; i++) 
      gxyz->SetPoint(i, eTrack.GetSegment(i)->X(), eTrack.GetSegment(i)->Y(), eTrack.GetSegment(i)->Z()); 
    gxyz->SetName("gxyz");
    gxyz->SetTitle("track: X vs Y vs Z");
    gxyz->SetMarkerStyle(21);
    gxyz->Draw("P");

    TGraphErrors *gtx = new TGraphErrors(nseg);
    for(int i=0; i<nseg; i++) 
      gtx->SetPoint(i, eTrack.GetSegment(i)->PID(), eTrack.GetSegment(i)->TX()); 
    TGraphErrors *gty = new TGraphErrors(nseg);
    for(int i=0; i<nseg; i++) 
      gty->SetPoint(i, eTrack.GetSegment(i)->PID(), eTrack.GetSegment(i)->TY()); 


    TVirtualPad *vp;
    vp = c1->cd(5);
    vp->SetGrid();
    gtx->SetLineColor(kBlue);
    gtx->SetMarkerStyle(24);
    gtx->Draw("ALP");

    vp = c1->cd(6);
    vp->SetGrid();
    gty->SetMarkerStyle(24);
    gty->SetLineColor(kRed);
    gty->Draw("ALP");
      }

}

//________________________________________________________________________________________
int EdbMomentumEstimator::PMSang_base(EdbTrackP &tr)
{
  // Version rewised by VT 13/05/2008
  //
  // Momentum estimation by multiple scattering (Annecy algorithm Oct-2007)
  //
  // Input: tr  - can be modified by the function
  //
  // calculate momentum in transverse and in longitudinal projections using the different 
  // measurements errors parametrisation
  // return value:  -99 - estimation impossible; 
  //                  0 - fit is not successful; 
  //                  1 - only one momentum component is fitted well
  //                  2 - both components are successfully fitted
  // "base" is stay for the original version - to be tested in comparison to the "final" version

  int nseg = tr.N();
  if(nseg<2)   { Log(1,"PMSang_base","Warning! nseg<2 (%d)- impossible estimate momentum!",nseg);             return -99;}
  int npl = tr.Npl();
  if(npl<nseg) { Log(1,"PMSang_base","Warning! npl<nseg (%d, %d) - use track.SetCounters() first",npl,nseg);  return -99;}
  int plmax = Max( tr.GetSegmentFirst()->PID(), tr.GetSegmentLast()->PID() ) + 1;
  if(plmax<1||plmax>1000)   { Log(1,"PMSang_base","Warning! plmax = %d - correct the segments PID's!",plmax); return -99;}
  Log(3,"PMSang_base","estimate track with %d segments %d plates",tr.N(), tr.Npl());

  float xmean,ymean,zmean,txmean,tymean,wmean;
  FitTrackLine(tr,xmean,ymean,zmean,txmean,tymean,wmean);    // calculate mean track parameters
  EdbAffine2D aff;
  aff.ShiftX(-xmean);
  aff.ShiftY(-ymean);
  aff.Rotate( -ATan2(tymean,txmean) );                       // rotate track to get longitudinal as tx, transverse as ty angle
  tr.Transform(aff);
  FitTrackLine(tr,xmean,ymean,zmean,txmean,tymean,wmean);    // calculate mean track parameters

  int minentr  = eMinEntr;        // min number of entries in the cell to accept the cell for fitting
  int stepmax  = npl/minentr;     // max step
  int size     = stepmax+1;       // vectors size

  TVectorF dax(size), day(size);
  TArrayI  nentr(size);

  Log(3,"PMSang_base","stepmax = %d",stepmax);

  EdbSegP *s1,*s2;
  for(int ist=1; ist<=stepmax; ist++)         // cycle by the step size
    {
      for(int i1=0; i1<nseg-1; i1++)          // cycle by the first seg
	{
	  s1 = tr.GetSegment(i1);
	  for(int i2=i1+1; i2<nseg; i2++)      // cycle by the second seg
	    {
	      s2 = tr.GetSegment(i2);
	      int icell = Abs(s2->PID()-s1->PID());
	      if( icell == ist ) {
		dax[icell-1]   += ( (ATan(s2->TX())- ATan(s1->TX())) * (ATan(s2->TX())- ATan(s1->TX())) );
		day[icell-1]   += ( (ATan(s2->TY())- ATan(s1->TY())) * (ATan(s2->TY())- ATan(s1->TY())) );
		nentr[icell-1] +=1;
	      }
	    }
	}
    }
 
  float maxX =0;                                  // maximum value for the function fit
  TVector vind(size), errvind(size);
  TVector errdax(size), errday(size);
  int ist=0;                                      // use the counter for case of missing cells 
  for(int i=0; i<size; i++) 
    {
      if( nentr[i] >= minentr ) {
	vind[ist]    = i+1;                           // x-coord is defined as the number of cells
	dax[ist]     = Sqrt( dax[ist]/nentr[i] );
	day[ist]     = Sqrt( day[ist]/nentr[i] );
	errvind[ist] = 0.25;
	errdax[ist]  = dax[ist]/CellWeight(npl,i+1);
	errday[ist]  = day[ist]/CellWeight(npl,i+1);
	maxX         = vind[ist];
	ist++;
      }
    }

  float dtx = GetDTx(txmean);  // measurements errors parametrization, longtudinal
  dtx*=dtx;
  float dty = GetDTy(tymean);  // measurements errors parametrization, transversal
  dty*=dty;

  float Zcorr = Sqrt(1+txmean*txmean+tymean*tymean);
  float x0    = eX0/1000/Zcorr;                       // the effective rad length in [mm]
 
  SafeDelete(eF1X);
  SafeDelete(eF1Y);
  SafeDelete(eGX);
  SafeDelete(eGY);
  
  const char *fitopt = "MQR"; //MQR

  bool statFitPX = false, statFitPY  = false;
  float initP = 1., minP=0., maxP=100.;                             // starting value for momentum in GeV

  eF1X = MCSErrorFunction_base("eF1X",x0,dtx);    eF1X->SetRange(0,maxX);
  eF1X->SetParameter(0, initP);
  //  eF1X->SetParameter(1, 0.002);
  eF1X->SetParLimits(0, minP, maxP);
  eGX=new TGraphErrors(vind,dax,errvind,errdax);
  eGX->Fit("eF1X",fitopt);
  ePx=eF1X->GetParameter(0);
  if( Abs(ePx-initP)<0.00001 ) {
    eF1X->SetParameter(0, 2*initP);
    eGX->Fit("eF1X",fitopt);
    ePx=eF1X->GetParameter(0);
    if( Abs(ePx - 2*initP)>0.00001 ) statFitPX=true;
  }
  else  statFitPX=true;

  if(statFitPX) eDPx=eF1X->GetParError(0);
  else { eDPx =-99; ePx = -99; }

  eF1Y = MCSErrorFunction_base("eF1Y",x0,dty);    eF1Y->SetRange(0,maxX);
  eF1Y->SetParameter(0,initP);
  eF1Y->SetParLimits(0, minP, maxP);
  eGY=new TGraphErrors(vind,day,errvind,errday);
  eGY->Fit("eF1Y",fitopt);
  ePy=eF1Y->GetParameter(0);
  if( Abs(ePy-initP)<0.00001 ) {
    eF1Y->SetParameter(0, 2*initP);
    eGY->Fit("eF1Y",fitopt);
    ePy=eF1Y->GetParameter(0);
    if( Abs(ePy - 2*initP)>0.00001 ) statFitPY=true;
  }
  else  statFitPY=true;
  if(statFitPY) eDPy=eF1Y->GetParError(0);
  else { eDPy =-99; ePy = -99; }

  EstimateMomentumError( ePx, npl, txmean, ePXmin, ePXmax );
  EstimateMomentumError( ePy, npl, tymean, ePYmin, ePYmax );

  float wx = statFitPX? 1./eDPx/eDPx : 0;
  float wy = statFitPY? 1./eDPy/eDPy : 0;
  if(statFitPX||statFitPY) {
    eP  = (ePx*wx + ePy*wy)/(wx+wy);   // TODO: check on MC the biases of different estimations
    eDP = Sqrt(1./(wx+wy));
  } else {eP = -1; eDP=-1; }

  if(gEDBDEBUGLEVEL>1)
    printf("id=%6d (%2d/%2d) px=%7.2f +-%5.2f (%6.2f : %6.2f)    py=%7.2f +-%5.2f  (%6.2f : %6.2f)  pmean =%7.2f  %d %d\n",
	   tr.ID(),npl,nseg,ePx,eDPx,ePXmin, ePXmax,ePy,eDPy,ePYmin, ePYmax, eP, statFitPX, statFitPY);
  
  return statFitPX+statFitPY;
}

//________________________________________________________________________________________
TF1 *EdbMomentumEstimator::MCSErrorFunction_base(const char *name, float x0, float dtx)
{
  //        dtx - the plane angle measurement error
  // return the function of the expected angular deviation vs range
  //
  // use the Highland-Lynch-Dahl formula for theta_rms_plane = 13.6 MeV/bcp*z*sqrt(x/x0)*(1+0.038*log(x/x0))  (PDG)
  // so the expected measured angle is sqrt( theta_rms_plane**2 + dtx**2)
  //
  // The constant term im the scattering formula is not 13.6 but 14.64, which
  // is the right reevaluated number, due to a calculation with the moliere
  // distribution. 13.6 is an approximation. See Geant3 or 4 references for more explanations.???????
  //
  // err(x) = sqrt(k*x*(1+0.038*log(x/x0))/p**2 + dtx)

  //  float k   = 14.64*14.64/x0;
  // 14.64*14.64/1000/1000 = 0.0002143296  - we need p in GeV
  // 13.6*13.6/1000/1000   = 0.0001849599  - we need p in GeV

  return new TF1(name,Form("sqrt(0.0002143296*x/%f*((1+0.038*log(x/(%f)))**2)/([0])**2+%f)",x0,x0,dtx));
  //return new TF1(name,Form("sqrt(0.0001849599*x/%f*(1+0.038*log(x/(%f)))/([0])**2+%f)",x0,x0,dtx));
}

//________________________________________________________________________________________
float EdbMomentumEstimator::P_MS(EdbTrackP &tr)
{
  // momentum estimation by multiple scattering (first test version (VT))

  int	 stepmax = 1;
  int    nms = 0;
  double tms = 0.;
  int ist = 0;

  float m;  // the mass of the particle
  eM<0? m = tr.M(): m=eM;

  EdbSegP *s1=0,*s2=0;

  double dx,dy,dz,ds;
  double dtx,dty,dts,fact,ax1,ax2,ay1,ay2,dax1,dax2,day1,day2;

  int nseg = tr.N(), i1 = 0, i2 = 0;

  for (ist=1; ist<=stepmax; ist++) {     // step size

    for (i1=0; i1<(nseg-ist); i1++) {       // for each step just once

      i2 = i1+ist;

      s1 = tr.GetSegment(i1);
      s2 = tr.GetSegment(i2);
	
      dx = s2->X()-s1->X();
      dy = s2->Y()-s1->Y();
      dz = s2->Z()-s1->Z();
      ds = Sqrt(dx*dx+dy*dy+dz*dz);
	
      ax1 = ATan(s1->TX());
      ax2 = ATan(s2->TX());
      ay1 = ATan(s1->TY());
      ay2 = ATan(s2->TY());
      dax1 = s1->STX();
      dax2 = s2->STX();
      day1 = s1->STY();
      day2 = s2->STY();
      dtx = (ax2-ax1);
      dty = (ay2-ay1);
      dts = dtx*dtx+dty*dty;
      fact = 1.+0.038*Log(ds/eX0);
      dts = (dts-dax1-dax2-day1-day2)/ds/fact/fact;
      //	if (dts < 0.) dts = 0.;
      tms += dts;
      nms++;
    }
  }

  if(tms<=0) { 
    printf("P_MS: BAD estimation for track %d: tms=%g  nms=%d\n",tr.ID(),tms,nms);
    return 10;   // with correct parameters setting this problem is usually happend for hard tracks >=10 GeV
  }
  double pbeta = 0., pbeta2 = 0.;
  pbeta = Sqrt((double)nms/tms/eX0)*0.01923;
  pbeta2 = pbeta*pbeta;
  double p = 0.5*(pbeta2 + Sqrt(pbeta2*pbeta2 + 4.*pbeta2*m*m));
  if (p <= 0.)
    p = 0.;
  else
    p = Sqrt(p);
  
  if (eDE_correction)
    {
      double dtot = 0., eTPb = 1000./1300., e = 0., tkin = 0.;
      s1 = tr.GetSegment(0);
      s2 = tr.GetSegment(nseg-1);

      dx = s2->X()-s1->X();
      dy = s2->Y()-s1->Y();
      dz = s2->Z()-s1->Z();
    
      dtot = Sqrt(dx*dx+dy*dy+dz*dz)*eTPb;

      double DE = EdbPhysics::DeAveragePb(p, m, dtot);
      tkin = Sqrt(p*p + m*m) - m;

      if (tkin < DE)
	{
	  tkin = 0.5*DE;
	  e = tkin + m;
	  p = Sqrt(e*e - m*m);
	  DE = EdbPhysics::DeAveragePb(p, m, dtot);
	}
      tkin = tkin + 0.5*DE;
      e = tkin + m;
      p = Sqrt(e*e - m*m);
    }
  
  return (float)p;
}

//____________________________________________________________________________________
int EdbMomentumEstimator::PMSang_base_A(EdbTrackP &tr)
{
  // Version revised by Andrea Russo 13/03/2009 based on PMSang_base() by VT
  //
  // include asymmetrical errors in scattering VS ncell graphs, based on ChiSquare distribution
  //
  // improved check on correctness of fit result, based on chisquare cut (FitProbability > 0.05)

  int nseg = tr.N();
  if(nseg<2)   { Log(1,"PMSang_base","Warning! nseg<2 (%d)- impossible estimate momentum!",nseg);             return -99;}
  int npl = tr.Npl();
  if(npl<nseg) { Log(1,"PMSang_base","Warning! npl<nseg (%d, %d) - use track.SetCounters() first",npl,nseg);  return -99;}
  int plmax = Max( tr.GetSegmentFirst()->PID(), tr.GetSegmentLast()->PID() ) + 1;
  if(plmax<1||plmax>1000)   { Log(1,"PMSang_base","Warning! plmax = %d - correct the segments PID's!",plmax); return -99;}
  Log(3,"PMSang_base","estimate track with %d segments %d plates",tr.N(), tr.Npl());

  float xmean,ymean,zmean,txmean,tymean,wmean;
  FitTrackLine(tr,xmean,ymean,zmean,txmean,tymean,wmean);    // calculate mean track parameters
  EdbAffine2D aff;
  aff.ShiftX(-xmean);
  aff.ShiftY(-ymean);
  aff.Rotate( -ATan2(tymean,txmean) );                       // rotate track to get longitudinal as tx, transverse as ty angle
  tr.Transform(aff);
  FitTrackLine(tr,xmean,ymean,zmean,txmean,tymean,wmean);    // calculate mean track parameters

  int minentr  = eMinEntr;        // min number of entries in the cell to accept the cell for fitting
  int stepmax  = npl/minentr;     // max step
  int size     = stepmax+1;       // vectors size

  TVectorF dax(size), day(size);
  TArrayI  nentr(size);

  Log(3,"PMSang_base_A","stepmax = %d",stepmax);


  EdbSegP *s1,*s2;
  for(int ist=1; ist<=stepmax; ist++)         // cycle by the step size
    {
      for(int i1=0; i1<nseg-1; i1++)          // cycle by the first seg
    	{
    	  s1 = tr.GetSegment(i1);
    	  for(int i2=i1+1; i2<nseg; i2++)      // cycle by the second seg
    	  {
    	      s2 = tr.GetSegment(i2);
    	      int icell = Abs(s2->PID()-s1->PID());
    	      if( icell == ist ) {
          		dax[icell-1]   += ( (ATan(s2->TX())- ATan(s1->TX())) * (ATan(s2->TX())- ATan(s1->TX())) );
          		day[icell-1]   += ( (ATan(s2->TY())- ATan(s1->TY())) * (ATan(s2->TY())- ATan(s1->TY())) );
          		nentr[icell-1] +=1;
    	      }
    	  }
    	}
    }

  float maxX =0;                                  // maximum value for the function fit
  TVector vind(size), errvind(size);
  
  TVector errdaxL(size), errdayL(size);           // L stands for low, lower error bar
  TVector errdaxH(size), errdayH(size);           // H stands for high, hogher error bar
  int ist=0;                                      // use the counter for case of missing cells 
  for(int i=0; i<size; i++) 
  {
    if( nentr[i] >= minentr ) {
        	float ndf = CellWeight(npl,i+1);           // CellWeight is interpreted as the ndf in the determination of dax and day

        	vind[ist]    = i+1;                           // x-coord is defined as the number of cells
        	dax[ist]     = Sqrt( dax[ist]/nentr[i] );
        	day[ist]     = Sqrt( day[ist]/nentr[i] );

        	errvind[ist] = 0.25;

        	//errdax[ist]  = dax[ist]/CellWeight(npl,i+1);
        	//errday[ist]  = day[ist]/CellWeight(npl,i+1);

        	errdaxL[ist]  = dax[ist] - dax[ist]/(Sqrt(TMath::ChisquareQuantile(0.84,ndf)/ndf));
        	errdaxH[ist]  = dax[ist]/(Sqrt(TMath::ChisquareQuantile(0.16,ndf)/ndf)) - dax[ist];

        	errdayL[ist]  = day[ist] - day[ist]/(Sqrt(TMath::ChisquareQuantile(0.84,CellWeight(npl,i+1))/ndf));
        	errdayH[ist]  = day[ist]/(Sqrt(TMath::ChisquareQuantile(0.16,ndf)/ndf)) - day[ist];

        	maxX         = vind[ist];
        	ist++;
    }
  }

  float dtx = GetDTx(txmean);  // measurements errors parametrization     
  dtx*=dtx;
  float dty = GetDTy(txmean);  // measurements errors parametrization
  dty*=dty;

  float Zcorr = Sqrt(1+txmean*txmean+tymean*tymean);
  float x0    = eX0/1000/Zcorr;                       // the effective rad length in [mm]

  SafeDelete(eF1X);
  SafeDelete(eF1Y);
  SafeDelete(eGAX);
  SafeDelete(eGAY);

  bool statFitPX = false, statFitPY  = false;
  float initP = 1., minP=0., maxP=100.;                             // starting value for momentum in GeV

  eF1X = MCSErrorFunction_base("eF1X",x0,dtx);    eF1X->SetRange(0,maxX);
  eF1X->SetParameter(0, initP);
  eF1X->SetParameter(1, 0.003);
  eF1X->SetParLimits(0, minP, maxP);
  
  eGAX=new TGraphAsymmErrors(vind,dax,errvind,errvind,errdaxL,errdaxH);
  
  const char *fitopt = "MQR"; //MQR
  eGAX->Fit("eF1X",fitopt);
  ePx=eF1X->GetParameter(0);
  
  /*
  if( Abs(ePx-initP)<0.00001 ) {
    eF1X->SetParameter(0, 2*initP);
    eGAX->Fit("eF1X",fitopt);
    ePx=eF1X->GetParameter(0);
    if( Abs(ePx - 2*initP)>0.00001 ) statFitPX=true;
  }
  else  statFitPX=true;
  */
  
  if(eF1X->GetChisquare() > TMath::ChisquareQuantile(0.05,eF1X->GetNDF())) //checking if fit procedure gives a reasonable result (i.e., p>0.05)
    statFitPX=false;
  else
    statFitPX=true;
  
  //statFitPX=true;


  if(statFitPX) eDPx=eF1X->GetParError(0);
  else { eDPx =-99; ePx = -99; }

  eF1Y = MCSErrorFunction_base("eF1Y",x0,dty);    eF1Y->SetRange(0,maxX);
  eF1Y->SetParameter(0,initP);
  eF1Y->SetParLimits(0, minP, maxP);

  eGAY=new TGraphAsymmErrors(vind,day,errvind,errvind,errdayL,errdayH);
  
  eGAY->Fit("eF1Y", fitopt);
  ePy=eF1Y->GetParameter(0);
  /*
  if( Abs(ePy-initP)<0.00001 ) {
    eF1Y->SetParameter(0, 2*initP);
    eGAY->Fit("eF1Y",fitopt);
    ePy=eF1Y->GetParameter(0);
    if( Abs(ePy - 2*initP)>0.00001 ) statFitPY=true;
  }
  else  statFitPY=true;
  */

  if(eF1Y->GetChisquare()> TMath::ChisquareQuantile(0.05,eF1Y->GetNDF())) //checking if fit procedure gives a reasonable result (i.e., p>0.05)
    statFitPY=false;
  else
    statFitPY=true;
  
  //statFitPY=true;

  if(statFitPY) eDPy=eF1Y->GetParError(0);
  else { eDPy =-99; ePy = -99; }

  EstimateMomentumError( ePx, npl, txmean, ePXmin, ePXmax );
  EstimateMomentumError( ePy, npl, tymean, ePYmin, ePYmax );

  float wx = statFitPX? 1./eDPx/eDPx : 0;
  float wy = statFitPY? 1./eDPy/eDPy : 0;
  if(statFitPX||statFitPY) {
    eP  = (ePx*wx + ePy*wy)/(wx+wy);   // TODO: check on MC the biases of different estimations
    eDP = Sqrt(1./(wx+wy));
  } else {eP = -1; eDP=-1; }

  if(gEDBDEBUGLEVEL>1)
    printf("id=%6d (%2d/%2d) px=%7.2f +-%5.2f (%6.2f : %6.2f)    py=%7.2f +-%5.2f  (%6.2f : %6.2f)  pmean =%7.2f  %d %d\n",
	   tr.ID(),npl,nseg,ePx,eDPx,ePXmin, ePXmax,ePy,eDPy,ePYmin, ePYmax, eP, statFitPX, statFitPY);
  
  
  return statFitPX+statFitPY;
}
