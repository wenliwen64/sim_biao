#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TMath.h>

using namespace std;
using namespace TMath;


Double_t cycle(Double_t Phi)
{
  if(Phi>2*Pi()) return Phi-2*Pi();
  if(Phi<0) return Phi+2*Pi();
  return Phi;
}

Double_t RoPhi(Double_t Phi)
{  
  TRandom3 flip(0);
  if(flip.Rndm()<0.5)
    {
      if(Phi>Pi())
	return Phi-Pi();
      else 
	return Phi+Pi();
    }
  else return Phi;
}

Int_t Phi_bin(Double_t Phi)
{
  if(Phi>=0&&Phi<0.5*TMath::Pi())
    {
      return 1;
    }
  else if(Phi>=0.5*TMath::Pi()&&Phi<TMath::Pi())
    {
      return 2;
    }
  else if(Phi>=TMath::Pi()&&Phi<1.5*TMath::Pi()) 
    {
      return 3;
    }
  else 
    {
      return 4;
    }
}


/*
Int_t CountNum(string ChWant,int charge,double Phi,double* Num,string Mode,double EPangle,double wedge)
{
  Phi=cycle(Phi-EPangle);
  if(Mode=="UD")
    {     
      if(ChWant=="+")
	{
	  if(charge>0)
	    {
	      if(fabs(Phi-0.5*Pi())<wedge) Num[0]+=charge;
	      if(fabs(Phi+0.5*Pi())<wedge) Num[1]+=charge;
	    }
	}
      if(ChWant=="-")
	{
	  if(charge<0)
	    {
	      if(fabs(Phi-0.5*Pi())<wedge) Num[0]+=charge;
	      if(fabs(Phi+0.5*Pi())<wedge) Num[1]+=charge;
	    }
	}
    }
  
  if(Mode=="LR")
    {     
      if(ChWant=="+")
	{
	  if(charge>0)
	    {
	      if((Phi-Pi())>-1*wedge||(Phi+Pi())<wedge) Num[0]+=charge;
	      if(fabs(Phi)<wedge) Num[1]+=charge;
	    }
	}
      if(ChWant=="-")
	{
	  if(charge<0)
	    {
	      if((Phi-Pi())>-1*wedge||(Phi+Pi())<wedge) Num[0]+=charge;
	      if(fabs(Phi)<wedge) Num[1]+=charge;
	    }
	}
    }
  return 0;
}
*/

Double_t CalA(Double_t* Num)
{
  if(Num[0]+Num[1]==0)
    {      
      return 9999;
    }
  else return (Num[0]-Num[1])/(Num[0]+Num[1]);
}


int main(int argc=0)
{
  
  Double_t v2=0.05;
  Double_t a=0.02;

  //***********************************************************

  char tmp[128];

  TH1D* hphi_p_east = new TH1D("phi_p_east","",100,-1,7); 
  TH1D* hphi_n_east = new TH1D("phi_n_east","",100,-1,7); 
  TH1D* hphi_p_west = new TH1D("phi_p_west","",100,-1,7); 
  TH1D* hphi_n_west = new TH1D("phi_n_west","",100,-1,7); 
 
  //*********************0 RP****1 EP************************

  TH1D* hpp_east[2];
  TH1D* hnn_east[2];
  TH1D* hpn_east[2];

  TH1D* hpp_west[2];
  TH1D* hnn_west[2]; 
  TH1D* hpn_west[2];

  TH1D* hpp[2];
  TH1D* hnn[2];
  TH1D* hpn[2];

  TH1D* hv2[2];

  for(Int_t i=0;i<2;i++)
    {
      sprintf(tmp,"pp_east_%d",i);
      hpp_east[i]= new TH1D(tmp,"",100,-0.1,0.1);
      sprintf(tmp,"nn_east_%d",i);
      hnn_east[i]= new TH1D(tmp,"",100,-0.1,0.1);
      sprintf(tmp,"pn_east_%d",i);
      hpn_east[i]= new TH1D(tmp,"",100,-0.1,0.1);


      sprintf(tmp,"pp_west_%d",i);
      hpp_west[i]= new TH1D(tmp,"",100,-0.1,0.1);
      sprintf(tmp,"nn_west_%d",i);
      hnn_west[i]= new TH1D(tmp,"",100,-0.1,0.1);
      sprintf(tmp,"pn_west_%d",i);
      hpn_west[i]= new TH1D(tmp,"",100,-0.1,0.1);


      sprintf(tmp,"pp_%d",i);
      hpp[i]= new TH1D(tmp,"",100,-0.1,0.1);
      sprintf(tmp,"nn_%d",i);
      hnn[i]= new TH1D(tmp,"",100,-0.1,0.1);
      sprintf(tmp,"pn_%d",i);
      hpn[i]= new TH1D(tmp,"",100,-0.1,0.1);

      sprintf(tmp,"v2_full_%d",i);
      hv2[i]= new TH1D(tmp,"",1000,-1,1);
    }

  
  TH1D* hRes = new TH1D("Res","",200,-1,1);

  TH1D* hEP[3];

  //******************0 east****1 west****2 full******************

  for(Int_t i=0;i<3;i++)
    {
      sprintf(tmp,"EP_%d",i);
      hEP[i]=new TH1D(tmp,"",1000,-1*TMath::Pi(),TMath::Pi());
    }

  TH1D* hv2ch[4];
  TH1D* hDelta[4];
  TH1D* hASqChPoUD[4]; 
  TH1D* hASqChNeUD[4];
  TH1D* hASqChPoLR[4]; 
  TH1D* hASqChNeLR[4];
  TH1D* hApAmChUD[4];
  TH1D* hApAmChLR[4];

  TH1D* hASqChPoStUD[4];
  TH1D* hASqChNeStUD[4];
  TH1D* hASqChPoStLR[4];
  TH1D* hASqChNeStLR[4];  
  TH1D* hApAmChStUD[4];
  TH1D* hApAmChStLR[4];

  TH2D* ChPlSqUD[4];
  TH2D* ChPlSqStUD[4];

  TH2D* ChMiSqUD[4];
  TH2D* ChMiSqStUD[4];
  
  TH2D* ChPlMiUD[4];
  TH2D* ChPlMiStUD[4];

  TH2D* ChPlSqLR[4];
  TH2D* ChPlSqStLR[4];

  TH2D* ChMiSqLR[4];
  TH2D* ChMiSqStLR[4];
  
  TH2D* ChPlMiLR[4];
  TH2D* ChPlMiStLR[4];


  TH1D* hanglediff_eastv2_westep = new TH1D("angle diff east v2 west ep", 100, -3.1415926, 3.1415926);
  TH1D* hanglediff_westv2_eastep = new TH1D("angle diff west v2 east ep", 100, -3.1415926, 3.1415926);
  //*******0 east_rp****1 west_rp****2 east_ep****3 west_ep*******

  for(Int_t i=0;i<4;i++)
    {
      
      sprintf(tmp,"v2ch_%d",i);
      hv2ch[i]=new TH1D(tmp,"",1000,-1.0,1.0);	  	  
      sprintf(tmp,"Delta_%d",i);
      hDelta[i]=new TH1D(tmp,"",2000,-1,1);
         
      sprintf(tmp,"ASqChPoUD_%d",i);
      hASqChPoUD[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ASqChNeUD_%d",i);
      hASqChNeUD[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ASqChPoLR_%d",i);
      hASqChPoLR[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ASqChNeLR_%d",i);
      hASqChNeLR[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ApAmChUD_%d",i);
      hApAmChUD[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ApAmChLR_%d",i);
      hApAmChLR[i]=new TH1D(tmp,"",2000,-1,1);
      
      sprintf(tmp,"ASqChPoStUD_%d",i);
      hASqChPoStUD[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ASqChNeStUD_%d",i);
      hASqChNeStUD[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ASqChPoStLR_%d",i);
      hASqChPoStLR[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ASqChNeStLR_%d",i);
      hASqChNeStLR[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ApAmChStUD_%d",i);
      hApAmChStUD[i]=new TH1D(tmp,"",2000,-1,1);
      sprintf(tmp,"ApAmChStLR_%d",i);
      hApAmChStLR[i]=new TH1D(tmp,"",2000,-1,1);
	  
      sprintf(tmp,"ChPlSqUD_%d",i);
      ChPlSqUD[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      sprintf(tmp,"ChPlSqStUD_%d",i);
      ChPlSqStUD[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      
      sprintf(tmp,"ChMiSqUD_%d",i);
      ChMiSqUD[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      sprintf(tmp,"ChMiSqStUD_%d",i);
      ChMiSqStUD[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      
      sprintf(tmp,"ChPlMiUD_%d",i);
      ChPlMiUD[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      sprintf(tmp,"ChPlMiStUD_%d",i);
      ChPlMiStUD[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      
      sprintf(tmp,"ChPlSqLR_%d",i);
      ChPlSqLR[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      sprintf(tmp,"ChPlSqStLR_%d",i);
      ChPlSqStLR[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
	  
      sprintf(tmp,"ChMiSqLR_%d",i);
      ChMiSqLR[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      sprintf(tmp,"ChMiSqStLR_%d",i);
      ChMiSqStLR[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      
      sprintf(tmp,"ChPlMiLR_%d",i);
      ChPlMiLR[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);
      sprintf(tmp,"ChPlMiStLR_%d",i);
      ChPlMiStLR[i]=new TH2D(tmp,"",200,-1,1,200,-1,1);    
    }  

  //************************************************************
  const Int_t nn=500;

  Int_t nEvent=0;

  //*******************************generate particles***********************************

  TRandom3 a1(0),a2(0),a3(0);
  Double_t b1,b2,b3;
  
  while(nEvent<1000)
    {
      Int_t nParticles=0;  

      Double_t Phi[nn];
      Int_t Charge[nn];
      Int_t Eta[nn];   

      Double_t EP_full=0;
      Double_t EP_east=0;
      Double_t EP_west=0; 
      Double_t EP_sin_east=0;
      Double_t EP_cos_east=0;
      Double_t EP_sin_west=0;
      Double_t EP_cos_west=0;

//===========Loop below is to generate particles and store eta information=========================
      while(nParticles<nn)
	{
	  //b1=a1.Rndm();
	  if(nParticles<0.5*nn)//positive
	    {	    
	      Int_t flag1=0;  
	      while(!flag1)
		{		 
		  b2=a2.Rndm();
		  b3=a3.Rndm();
		  if(b2<(2.0/3.0)*(1+2.0*v2*cos(4*TMath::Pi()*b3)+2.0*a*sin(2*TMath::Pi()*b3)))
		    {	    
		      Phi[nParticles]=2.0*TMath::Pi()*b3;
		      Charge[nParticles]=1;
		      flag1=1;		     
		      if(nParticles%2==0) //=======So gurantee half eta > 0 half eta < 0==================
			{
			  Eta[nParticles]=-1;
			  hphi_p_east->Fill(Phi[nParticles]);
			  EP_sin_east+=sin(2*Phi[nParticles]);
			  EP_cos_east+=cos(2*Phi[nParticles]);
			}
		      else
			{
			  Eta[nParticles]=1;
			  hphi_p_west->Fill(Phi[nParticles]);
			  EP_sin_west+=sin(2*Phi[nParticles]);
			  EP_cos_west+=cos(2*Phi[nParticles]);
			}
		    }
		  else
		    {
		      flag1=0;
		    }
		}
	      nParticles++;
	    }	      
	 
	  else if(nParticles<nn)//negative
	    {
	      Int_t flag2=0;  
	      while(!flag2)
		{
		  b2=a2.Rndm();
		  b3=a3.Rndm();
		  if(b2<(2.0/3.0)*(1+2.0*v2*cos(4*TMath::Pi()*b3)-2.0*a*sin(2*TMath::Pi()*b3)))
		    {	    
		      Phi[nParticles]=2.0*TMath::Pi()*b3;
		      Charge[nParticles]=-1;
		      flag2=1;		      
		      if(nParticles%2==1)
			{
			  Eta[nParticles]=-1;
			  hphi_n_east->Fill(Phi[nParticles]);
			  EP_sin_east+=sin(2*Phi[nParticles]);
			  EP_cos_east+=cos(2*Phi[nParticles]);
			}
		      else 
			{
			  Eta[nParticles]=1;
			  hphi_n_west->Fill(Phi[nParticles]);
			  EP_sin_west+=sin(2*Phi[nParticles]);
			  EP_cos_west+=cos(2*Phi[nParticles]);
			}
		    }
		  else
		    {
		      flag2=0;
		    }
		}
	      nParticles++;
	    }
	}

     

      //******************************Event Plane*******************************

      EP_east=atan2(EP_sin_east,EP_cos_east)/2.0;
      EP_west=atan2(EP_sin_west,EP_cos_west)/2.0;
      EP_full=atan2(EP_sin_west+EP_sin_east,EP_cos_west+EP_cos_east)/2.0;

      hEP[0]->Fill(EP_east);
      hEP[1]->Fill(EP_west);
      hEP[2]->Fill(EP_full);
   

      hRes->Fill(cos(2*(EP_east-EP_west)));
   
      //************************************************************************
      
      //***************************A and gamma**********************************
     
      Int_t p_num_east=0,p_num_west=0;
      Int_t n_num_east=0,n_num_west=0;


      //**********************RP************************

      double *NumOfChPoUD_east_rp=new double[2]; //======Number of positively charged particles up-0 and down-1 in east sub-event========
      double *NumOfChPoStUD_east_rp=new double[2]; //=======St??=======
      double *NumOfChNeUD_east_rp=new double[2];
      double *NumOfChNeStUD_east_rp=new double[2];
      
      double *NumOfChPoLR_east_rp=new double[2];
      double *NumOfChPoStLR_east_rp=new double[2];
      double *NumOfChNeLR_east_rp=new double[2];
      double *NumOfChNeStLR_east_rp=new double[2];
      
      double *NumOfChPoUD_west_rp=new double[2];
      double *NumOfChPoStUD_west_rp=new double[2];
      double *NumOfChNeUD_west_rp=new double[2];
      double *NumOfChNeStUD_west_rp=new double[2];
      
      double *NumOfChPoLR_west_rp=new double[2];
      double *NumOfChPoStLR_west_rp=new double[2];
      double *NumOfChNeLR_west_rp=new double[2];
      double *NumOfChNeStLR_west_rp=new double[2];

      double p_sumcos_east_rp,p_sumcos_west_rp; //======sum of positive particle cosine terms======
      double p_sumsin_east_rp,p_sumsin_west_rp;
      double n_sumcos_east_rp,n_sumcos_west_rp;
      double n_sumsin_east_rp,n_sumsin_west_rp;
      
      double p_sumcos2_east_rp,p_sumcos2_west_rp;
      double p_sumsin2_east_rp,p_sumsin2_west_rp;
      double n_sumcos2_east_rp,n_sumcos2_west_rp;
      double n_sumsin2_east_rp,n_sumsin2_west_rp;

      Double_t p_v2sum_east_rp=0,p_v2sum_west_rp=0;
      Double_t n_v2sum_east_rp=0,n_v2sum_west_rp=0;      

      p_sumcos_east_rp=0;p_sumcos_west_rp=0;
      p_sumsin_east_rp=0;p_sumsin_west_rp=0;
      n_sumcos_east_rp=0;n_sumcos_west_rp=0;
      n_sumsin_east_rp=0;n_sumsin_west_rp=0;
      
      p_sumcos2_east_rp=0;p_sumcos2_west_rp=0;
      p_sumsin2_east_rp=0;p_sumsin2_west_rp=0;
      n_sumcos2_east_rp=0;n_sumcos2_west_rp=0;
      n_sumsin2_east_rp=0;n_sumsin2_west_rp=0;


      NumOfChPoUD_east_rp[0]=NumOfChPoUD_east_rp[1]=0;
      NumOfChPoStUD_east_rp[0]=NumOfChPoStUD_east_rp[1]=0;
      NumOfChNeUD_east_rp[0]=NumOfChNeUD_east_rp[1]=0;
      NumOfChNeStUD_east_rp[0]=NumOfChNeStUD_east_rp[1]=0;
      
      NumOfChPoLR_east_rp[0]=NumOfChPoLR_east_rp[1]=0;
      NumOfChPoStLR_east_rp[0]=NumOfChPoStLR_east_rp[1]=0;
      NumOfChNeLR_east_rp[0]=NumOfChNeLR_east_rp[1]=0;
      NumOfChNeStLR_east_rp[0]=NumOfChNeStLR_east_rp[1]=0;
      
      NumOfChPoUD_west_rp[0]=NumOfChPoUD_west_rp[1]=0;
      NumOfChPoStUD_west_rp[0]=NumOfChPoStUD_west_rp[1]=0;
      NumOfChNeUD_west_rp[0]=NumOfChNeUD_west_rp[1]=0;
      NumOfChNeStUD_west_rp[0]=NumOfChNeStUD_west_rp[1]=0;
      
      NumOfChPoLR_west_rp[0]=NumOfChPoLR_west_rp[1]=0;
      NumOfChPoStLR_west_rp[0]=NumOfChPoStLR_west_rp[1]=0;
      NumOfChNeLR_west_rp[0]=NumOfChNeLR_west_rp[1]=0;
      NumOfChNeStLR_west_rp[0]=NumOfChNeStLR_west_rp[1]=0;

    
      //**********************EP************************
  
      double *NumOfChPoUD_east_ep=new double[2];
      double *NumOfChPoStUD_east_ep=new double[2];
      double *NumOfChNeUD_east_ep=new double[2];
      double *NumOfChNeStUD_east_ep=new double[2];
      
      double *NumOfChPoLR_east_ep=new double[2];
      double *NumOfChPoStLR_east_ep=new double[2];
      double *NumOfChNeLR_east_ep=new double[2];
      double *NumOfChNeStLR_east_ep=new double[2];
      
      double *NumOfChPoUD_west_ep=new double[2];
      double *NumOfChPoStUD_west_ep=new double[2];
      double *NumOfChNeUD_west_ep=new double[2];
      double *NumOfChNeStUD_west_ep=new double[2];
      
      double *NumOfChPoLR_west_ep=new double[2];
      double *NumOfChPoStLR_west_ep=new double[2];
      double *NumOfChNeLR_west_ep=new double[2];
      double *NumOfChNeStLR_west_ep=new double[2];
      
      double p_sumcos_east_ep,p_sumcos_west_ep;
      double p_sumsin_east_ep,p_sumsin_west_ep;
      double n_sumcos_east_ep,n_sumcos_west_ep;
      double n_sumsin_east_ep,n_sumsin_west_ep;
      
      double p_sumcos2_east_ep,p_sumcos2_west_ep;
      double p_sumsin2_east_ep,p_sumsin2_west_ep;
      double n_sumcos2_east_ep,n_sumcos2_west_ep;
      double n_sumsin2_east_ep,n_sumsin2_west_ep;

      Double_t p_v2sum_east_ep=0,p_v2sum_west_ep=0;
      Double_t n_v2sum_east_ep=0,n_v2sum_west_ep=0;
      
      p_sumcos_east_ep=0;p_sumcos_west_ep=0;
      p_sumsin_east_ep=0;p_sumsin_west_ep=0;
      n_sumcos_east_ep=0;n_sumcos_west_ep=0;
      n_sumsin_east_ep=0;n_sumsin_west_ep=0;
      
      p_sumcos2_east_ep=0;p_sumcos2_west_ep=0;
      p_sumsin2_east_ep=0;p_sumsin2_west_ep=0;
      n_sumcos2_east_ep=0;n_sumcos2_west_ep=0;
      n_sumsin2_east_ep=0;n_sumsin2_west_ep=0;


      NumOfChPoUD_east_ep[0]=NumOfChPoUD_east_ep[1]=0;
      NumOfChPoStUD_east_ep[0]=NumOfChPoStUD_east_ep[1]=0;
      NumOfChNeUD_east_ep[0]=NumOfChNeUD_east_ep[1]=0;
      NumOfChNeStUD_east_ep[0]=NumOfChNeStUD_east_ep[1]=0;
      
      NumOfChPoLR_east_ep[0]=NumOfChPoLR_east_ep[1]=0;
      NumOfChPoStLR_east_ep[0]=NumOfChPoStLR_east_ep[1]=0;
      NumOfChNeLR_east_ep[0]=NumOfChNeLR_east_ep[1]=0;
      NumOfChNeStLR_east_ep[0]=NumOfChNeStLR_east_ep[1]=0;
      
      NumOfChPoUD_west_ep[0]=NumOfChPoUD_west_ep[1]=0;
      NumOfChPoStUD_west_ep[0]=NumOfChPoStUD_west_ep[1]=0;
      NumOfChNeUD_west_ep[0]=NumOfChNeUD_west_ep[1]=0;
      NumOfChNeStUD_west_ep[0]=NumOfChNeStUD_west_ep[1]=0;
      
      NumOfChPoLR_west_ep[0]=NumOfChPoLR_west_ep[1]=0;
      NumOfChPoStLR_west_ep[0]=NumOfChPoStLR_west_ep[1]=0;
      NumOfChNeLR_west_ep[0]=NumOfChNeLR_west_ep[1]=0;
      NumOfChNeStLR_west_ep[0]=NumOfChNeStLR_west_ep[1]=0;
    
     

      //***********************Full EP*****************************

      double p_sumcos_ep=0;
      double p_sumsin_ep=0;
      double n_sumcos_ep=0;
      double n_sumsin_ep=0;
      
      double p_sumcos2_ep=0;
      double p_sumsin2_ep=0;
      double n_sumcos2_ep=0;
      double n_sumsin2_ep=0;

      Double_t p_v2sum_ep=0;
      Double_t n_v2sum_ep=0;

      for(int i=0;i<nn;i++)
	{
	  
	  if(Eta[i]==-1)
	    {	     
	      if(Charge[i]==1)
		{
		  p_num_east++;

		  //**********************RP******************************

		  switch (Phi_bin(Phi[i]))
		    {
		    case 1: 
		      NumOfChPoUD_east_rp[0]++;NumOfChPoLR_east_rp[1]++;break;
		    case 2: 
		      NumOfChPoUD_east_rp[0]++;NumOfChPoLR_east_rp[0]++;break;
		    case 3:
		      NumOfChPoUD_east_rp[1]++;NumOfChPoLR_east_rp[0]++;break;
		    case 4:
		      NumOfChPoUD_east_rp[1]++;NumOfChPoLR_east_rp[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(Phi[i])))
		    {
		    case 1: 
		      NumOfChPoStUD_east_rp[0]++;NumOfChPoStLR_east_rp[1]++;break;
		    case 2: 
		      NumOfChPoStUD_east_rp[0]++;NumOfChPoStLR_east_rp[0]++;break;
		    case 3:
		      NumOfChPoStUD_east_rp[1]++;NumOfChPoStLR_east_rp[0]++;break;
		    case 4:
		      NumOfChPoStUD_east_rp[1]++;NumOfChPoStLR_east_rp[1]++;break;
		    }

		  p_sumcos_east_rp+=cos(Phi[i]);
		  p_sumsin_east_rp+=sin(Phi[i]);
		  p_sumcos2_east_rp+=cos(Phi[i])*cos(Phi[i]);
		  p_sumsin2_east_rp+=sin(Phi[i])*sin(Phi[i]);
		  p_v2sum_east_rp+=cos(2*Phi[i]);
			  
		  //**********************EP******************************

		  switch (Phi_bin(cycle(Phi[i]-EP_west)))
		    {
		    case 1: 
		      NumOfChPoUD_east_ep[0]++;NumOfChPoLR_east_ep[1]++;break;
		    case 2: 
		      NumOfChPoUD_east_ep[0]++;NumOfChPoLR_east_ep[0]++;break;
		    case 3:
		      NumOfChPoUD_east_ep[1]++;NumOfChPoLR_east_ep[0]++;break;
		    case 4:
		      NumOfChPoUD_east_ep[1]++;NumOfChPoLR_east_ep[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(cycle(Phi[i]-EP_west))))
		    {
		    case 1: 
		      NumOfChPoStUD_east_ep[0]++;NumOfChPoStLR_east_ep[1]++;break;
		    case 2: 
		      NumOfChPoStUD_east_ep[0]++;NumOfChPoStLR_east_ep[0]++;break;
		    case 3:
		      NumOfChPoStUD_east_ep[1]++;NumOfChPoStLR_east_ep[0]++;break;
		    case 4:
		      NumOfChPoStUD_east_ep[1]++;NumOfChPoStLR_east_ep[1]++;break;
		    }

		  p_sumcos_east_ep+=cos(Phi[i]-EP_west);
		  p_sumsin_east_ep+=sin(Phi[i]-EP_west);
		  p_sumcos2_east_ep+=cos(Phi[i]-EP_west)*cos(Phi[i]-EP_west);
		  p_sumsin2_east_ep+=sin(Phi[i]-EP_west)*sin(Phi[i]-EP_west);
		  p_v2sum_east_ep+=cos(2*(Phi[i]-EP_west));

		  //*************************Full EP*****************************

		  p_sumcos_ep+=cos(Phi[i]-EP_full);
		  p_sumsin_ep+=sin(Phi[i]-EP_full);
		  p_sumcos2_ep+=cos(Phi[i]-EP_full)*cos(Phi[i]-EP_full);
		  p_sumsin2_ep+=sin(Phi[i]-EP_full)*sin(Phi[i]-EP_full);
		  p_v2sum_ep+=cos(2*(Phi[i]-EP_full));	       		
		}
	      else 
		{
		  n_num_east++;

		  //**********************RP******************************

		  switch (Phi_bin(Phi[i]))
		    {
		    case 1: 
		      NumOfChNeUD_east_rp[0]++;NumOfChNeLR_east_rp[1]++;break;
		    case 2: 
		      NumOfChNeUD_east_rp[0]++;NumOfChNeLR_east_rp[0]++;break;
		    case 3:
		      NumOfChNeUD_east_rp[1]++;NumOfChNeLR_east_rp[0]++;break;
		    case 4:
		      NumOfChNeUD_east_rp[1]++;NumOfChNeLR_east_rp[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(Phi[i])))
		    {
		    case 1: 
		      NumOfChNeStUD_east_rp[0]++;NumOfChNeStLR_east_rp[1]++;break;
		    case 2: 
		      NumOfChNeStUD_east_rp[0]++;NumOfChNeStLR_east_rp[0]++;break;
		    case 3:
		      NumOfChNeStUD_east_rp[1]++;NumOfChNeStLR_east_rp[0]++;break;
		    case 4:
		      NumOfChNeStUD_east_rp[1]++;NumOfChNeStLR_east_rp[1]++;break;
		    }

		  n_sumcos_east_rp+=cos(Phi[i]);
		  n_sumsin_east_rp+=sin(Phi[i]);
		  n_sumcos2_east_rp+=cos(Phi[i])*cos(Phi[i]);
		  n_sumsin2_east_rp+=sin(Phi[i])*sin(Phi[i]);
		  n_v2sum_east_rp+=cos(2*Phi[i]);		

		  //**********************EP******************************

		  switch (Phi_bin(cycle(Phi[i]-EP_west)))
		    {
		    case 1: 
		      NumOfChNeUD_east_ep[0]++;NumOfChNeLR_east_ep[1]++;break;
		    case 2: 
		      NumOfChNeUD_east_ep[0]++;NumOfChNeLR_east_ep[0]++;break;
		    case 3:
		      NumOfChNeUD_east_ep[1]++;NumOfChNeLR_east_ep[0]++;break;
		    case 4:
		      NumOfChNeUD_east_ep[1]++;NumOfChNeLR_east_ep[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(cycle(Phi[i]-EP_west))))
		    {
		    case 1: 
		      NumOfChNeStUD_east_ep[0]++;NumOfChNeStLR_east_ep[1]++;break;
		    case 2: 
		      NumOfChNeStUD_east_ep[0]++;NumOfChNeStLR_east_ep[0]++;break;
		    case 3:
		      NumOfChNeStUD_east_ep[1]++;NumOfChNeStLR_east_ep[0]++;break;
		    case 4:
		      NumOfChNeStUD_east_ep[1]++;NumOfChNeStLR_east_ep[1]++;break;
		    }

		  n_sumcos_east_ep+=cos(Phi[i]-EP_west);
		  n_sumsin_east_ep+=sin(Phi[i]-EP_west);
		  n_sumcos2_east_ep+=cos(Phi[i]-EP_west)*cos(Phi[i]-EP_west);
		  n_sumsin2_east_ep+=sin(Phi[i]-EP_west)*sin(Phi[i]-EP_west);
		  n_v2sum_east_ep+=cos(2*(Phi[i]-EP_west));

		  //*************************Full EP*****************************

		  n_sumcos_ep+=cos(Phi[i]-EP_full);
		  n_sumsin_ep+=sin(Phi[i]-EP_full);
		  n_sumcos2_ep+=cos(Phi[i]-EP_full)*cos(Phi[i]-EP_full);
		  n_sumsin2_ep+=sin(Phi[i]-EP_full)*sin(Phi[i]-EP_full);
		  n_v2sum_ep+=cos(2*(Phi[i]-EP_full));		 
		}
	    }

	  else
	    {	     	     	
	      if(Charge[i]==1)
		{
		  p_num_west++;

		  //**********************RP******************************

		  switch (Phi_bin(Phi[i]))
		    {
		    case 1: 
		      NumOfChPoUD_west_rp[0]++;NumOfChPoLR_west_rp[1]++;break;
		    case 2: 
		      NumOfChPoUD_west_rp[0]++;NumOfChPoLR_west_rp[0]++;break;
		    case 3:
		      NumOfChPoUD_west_rp[1]++;NumOfChPoLR_west_rp[0]++;break;
		    case 4:
		      NumOfChPoUD_west_rp[1]++;NumOfChPoLR_west_rp[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(Phi[i])))
		    {
		    case 1: 
		      NumOfChPoStUD_west_rp[0]++;NumOfChPoStLR_west_rp[1]++;break;
		    case 2: 
		      NumOfChPoStUD_west_rp[0]++;NumOfChPoStLR_west_rp[0]++;break;
		    case 3:
		      NumOfChPoStUD_west_rp[1]++;NumOfChPoStLR_west_rp[0]++;break;
		    case 4:
		      NumOfChPoStUD_west_rp[1]++;NumOfChPoStLR_west_rp[1]++;break;
		    }
		  
		  p_sumcos_west_rp+=cos(Phi[i]);		 
		  p_sumsin_west_rp+=sin(Phi[i]);
		  p_sumcos2_west_rp+=cos(Phi[i])*cos(Phi[i]);
		  p_sumsin2_west_rp+=sin(Phi[i])*sin(Phi[i]);
		  p_v2sum_west_rp+=cos(2*Phi[i]);		 

		  //**********************EP******************************

		  switch (Phi_bin(cycle(Phi[i]-EP_east)))
		    {
		    case 1: 
		      NumOfChPoUD_west_ep[0]++;NumOfChPoLR_west_ep[1]++;break;
		    case 2: 
		      NumOfChPoUD_west_ep[0]++;NumOfChPoLR_west_ep[0]++;break;
		    case 3:
		      NumOfChPoUD_west_ep[1]++;NumOfChPoLR_west_ep[0]++;break;
		    case 4:
		      NumOfChPoUD_west_ep[1]++;NumOfChPoLR_west_ep[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(cycle(Phi[i]-EP_east))))
		    {
		    case 1: 
		      NumOfChPoStUD_west_ep[0]++;NumOfChPoStLR_west_ep[1]++;break;
		    case 2: 
		      NumOfChPoStUD_west_ep[0]++;NumOfChPoStLR_west_ep[0]++;break;
		    case 3:
		      NumOfChPoStUD_west_ep[1]++;NumOfChPoStLR_west_ep[0]++;break;
		    case 4:
		      NumOfChPoStUD_west_ep[1]++;NumOfChPoStLR_west_ep[1]++;break;
		    }
		  
		  p_sumcos_west_ep+=cos(Phi[i]-EP_east);		 
		  p_sumsin_west_ep+=sin(Phi[i]-EP_east);
		  p_sumcos2_west_ep+=cos(Phi[i]-EP_east)*cos(Phi[i]-EP_east);
		  p_sumsin2_west_ep+=sin(Phi[i]-EP_east)*sin(Phi[i]-EP_east);
		  p_v2sum_west_ep+=cos(2*(Phi[i]-EP_east));

		  //*************************Full EP*****************************

		  p_sumcos_ep+=cos(Phi[i]-EP_full);
		  p_sumsin_ep+=sin(Phi[i]-EP_full);
		  p_sumcos2_ep+=cos(Phi[i]-EP_full)*cos(Phi[i]-EP_full);
		  p_sumsin2_ep+=sin(Phi[i]-EP_full)*sin(Phi[i]-EP_full);
		  p_v2sum_ep+=cos(2*(Phi[i]-EP_full));
		}
	      else 
		{

		  n_num_west++;

		  //**********************RP******************************

		  switch (Phi_bin(Phi[i]))
		    {
		    case 1: 
		      NumOfChNeUD_west_rp[0]++;NumOfChNeLR_west_rp[1]++;break;
		    case 2: 
		      NumOfChNeUD_west_rp[0]++;NumOfChNeLR_west_rp[0]++;break;
		    case 3:
		      NumOfChNeUD_west_rp[1]++;NumOfChNeLR_west_rp[0]++;break;
		    case 4:
		      NumOfChNeUD_west_rp[1]++;NumOfChNeLR_west_rp[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(Phi[i])))
		    {
		    case 1: 
		      NumOfChNeStUD_west_rp[0]++;NumOfChNeStLR_west_rp[1]++;break;
		    case 2: 
		      NumOfChNeStUD_west_rp[0]++;NumOfChNeStLR_west_rp[0]++;break;
		    case 3:
		      NumOfChNeStUD_west_rp[1]++;NumOfChNeStLR_west_rp[0]++;break;
		    case 4:
		      NumOfChNeStUD_west_rp[1]++;NumOfChNeStLR_west_rp[1]++;break;
		    }

		  n_sumcos_west_rp+=cos(Phi[i]);
		  n_sumsin_west_rp+=sin(Phi[i]);
		  n_sumcos2_west_rp+=cos(Phi[i])*cos(Phi[i]);
		  n_sumsin2_west_rp+=sin(Phi[i])*sin(Phi[i]);
		  n_v2sum_west_rp+=cos(2*Phi[i]);
		 
		  //**********************EP******************************

		  switch (Phi_bin(cycle(Phi[i]-EP_east)))
		    {
		    case 1: 
		      NumOfChNeUD_west_ep[0]++;NumOfChNeLR_west_ep[1]++;break;
		    case 2: 
		      NumOfChNeUD_west_ep[0]++;NumOfChNeLR_west_ep[0]++;break;
		    case 3:
		      NumOfChNeUD_west_ep[1]++;NumOfChNeLR_west_ep[0]++;break;
		    case 4:
		      NumOfChNeUD_west_ep[1]++;NumOfChNeLR_west_ep[1]++;break;
		    }

		  switch (Phi_bin(RoPhi(cycle(Phi[i]-EP_east))))
		    {
		    case 1: 
		      NumOfChNeStUD_west_ep[0]++;NumOfChNeStLR_west_ep[1]++;break;
		    case 2: 
		      NumOfChNeStUD_west_ep[0]++;NumOfChNeStLR_west_ep[0]++;break;
		    case 3:
		      NumOfChNeStUD_west_ep[1]++;NumOfChNeStLR_west_ep[0]++;break;
		    case 4:
		      NumOfChNeStUD_west_ep[1]++;NumOfChNeStLR_west_ep[1]++;break;
		    }

		  n_sumcos_west_ep+=cos(Phi[i]-EP_east);
		  n_sumsin_west_ep+=sin(Phi[i]-EP_east);
		  n_sumcos2_west_ep+=cos(Phi[i]-EP_east)*cos(Phi[i]-EP_east);
		  n_sumsin2_west_ep+=sin(Phi[i]-EP_east)*sin(Phi[i]-EP_east);
		  n_v2sum_west_ep+=cos(2*(Phi[i]-EP_east));

		  //*************************Full EP*****************************

		  n_sumcos_ep+=cos(Phi[i]-EP_full);
		  n_sumsin_ep+=sin(Phi[i]-EP_full);
		  n_sumcos2_ep+=cos(Phi[i]-EP_full)*cos(Phi[i]-EP_full);
		  n_sumsin2_ep+=sin(Phi[i]-EP_full)*sin(Phi[i]-EP_full);
		  n_v2sum_ep+=cos(2*(Phi[i]-EP_full));		
		} 
	    }
	}


    
      //************************************************************************
      
      //****************************calculation*********************************

      //**********************************RP************************************

      double AChPoUD_east_rp,AChNeUD_east_rp;
      double ASqChPoUD_east_rp,ASqChNeUD_east_rp;
      double ApAmChUD_east_rp;
      
      double AChPoLR_east_rp,AChNeLR_east_rp;
      double ASqChPoLR_east_rp,ASqChNeLR_east_rp;
      double ApAmChLR_east_rp;
  
      double AChPoUD_west_rp,AChNeUD_west_rp;
      double ASqChPoUD_west_rp,ASqChNeUD_west_rp;
      double ApAmChUD_west_rp;
      
      double AChPoLR_west_rp,AChNeLR_west_rp;
      double ASqChPoLR_west_rp,ASqChNeLR_west_rp;
      double ApAmChLR_west_rp;

      double AChPoStUD_east_rp,AChNeStUD_east_rp;
      double ASqChPoStUD_east_rp,ASqChNeStUD_east_rp;
      double ApAmChStUD_east_rp;
      
      double AChPoStLR_east_rp,AChNeStLR_east_rp;
      double ASqChPoStLR_east_rp,ASqChNeStLR_east_rp;
      double ApAmChStLR_east_rp;
  
      double AChPoStUD_west_rp,AChNeStUD_west_rp;
      double ASqChPoStUD_west_rp,ASqChNeStUD_west_rp;
      double ApAmChStUD_west_rp;
      
      double AChPoStLR_west_rp,AChNeStLR_west_rp;
      double ASqChPoStLR_west_rp,ASqChNeStLR_west_rp;
      double ApAmChStLR_west_rp;   


      Double_t v2ch_east_rp,v2ch_west_rp,v2full_rp;

      v2ch_east_rp=(p_v2sum_east_rp+n_v2sum_east_rp)/(p_num_east+n_num_east);
      v2ch_west_rp=(p_v2sum_west_rp+n_v2sum_west_rp)/(p_num_west+n_num_west);
      v2full_rp=(p_v2sum_east_rp+n_v2sum_east_rp+p_v2sum_west_rp+n_v2sum_west_rp)/(p_num_east+n_num_east+p_num_west+n_num_west);   

      hv2ch[0]->Fill(v2ch_east_rp);
      hv2ch[1]->Fill(v2ch_west_rp);
      hv2[0]->Fill(v2full_rp);

      int iflag1=1;
      
//==============================????????========================
      if((AChPoUD_east_rp=CalA(NumOfChPoUD_east_rp))>0.8) {iflag1=0;}
      if((AChPoStUD_east_rp=CalA(NumOfChPoStUD_east_rp))>0.8) {iflag1=0;}
      if((AChNeUD_east_rp=CalA(NumOfChNeUD_east_rp))>0.8) {iflag1=0;}
      if((AChNeStUD_east_rp=CalA(NumOfChNeStUD_east_rp))>0.8) {iflag1=0;}
	 
      if((AChPoLR_east_rp=CalA(NumOfChPoLR_east_rp))>0.8) {iflag1=0;}
      if((AChPoStLR_east_rp=CalA(NumOfChPoStLR_east_rp))>0.8) {iflag1=0;}
      if((AChNeLR_east_rp=CalA(NumOfChNeLR_east_rp))>0.8) {iflag1=0;}
      if((AChNeStLR_east_rp=CalA(NumOfChNeStLR_east_rp))>0.8) {iflag1=0;}


      if(iflag1==1)
	{
	  ASqChPoUD_east_rp=pow(AChPoUD_east_rp,2);ASqChPoStUD_east_rp=pow(AChPoStUD_east_rp,2);
	  ASqChNeUD_east_rp=pow(AChNeUD_east_rp,2);ASqChNeStUD_east_rp=pow(AChNeStUD_east_rp,2);
	  ApAmChUD_east_rp=AChPoUD_east_rp*AChNeUD_east_rp;ApAmChStUD_east_rp=AChPoStUD_east_rp*AChNeStUD_east_rp;
	  	  
	  ASqChPoLR_east_rp=pow(AChPoLR_east_rp,2);ASqChPoStLR_east_rp=pow(AChPoStLR_east_rp,2);
	  ASqChNeLR_east_rp=pow(AChNeLR_east_rp,2);ASqChNeStLR_east_rp=pow(AChNeStLR_east_rp,2);
	  ApAmChLR_east_rp=AChPoLR_east_rp*AChNeLR_east_rp;ApAmChStLR_east_rp=AChPoStLR_east_rp*AChNeStLR_east_rp;
	     
	  Double_t delta_east_rp=(ASqChPoUD_east_rp-ASqChPoStUD_east_rp+ASqChNeUD_east_rp-ASqChNeStUD_east_rp)/2-(ASqChPoLR_east_rp-ASqChPoStLR_east_rp+ASqChNeLR_east_rp-ASqChNeStLR_east_rp)/2-(ApAmChUD_east_rp-ApAmChStUD_east_rp-ApAmChLR_east_rp+ApAmChStLR_east_rp);

	  hDelta[0]->Fill(delta_east_rp);	     

	  hASqChPoUD[0]->Fill(ASqChPoUD_east_rp);
	  hASqChNeUD[0]->Fill(ASqChNeUD_east_rp);
	  hASqChPoLR[0]->Fill(ASqChPoLR_east_rp);
	  hASqChNeLR[0]->Fill(ASqChNeLR_east_rp);
	  hApAmChUD[0]->Fill(ApAmChUD_east_rp);
	  hApAmChLR[0]->Fill(ApAmChLR_east_rp);
	  
	  hASqChPoStUD[0]->Fill(ASqChPoStUD_east_rp);
	  hASqChNeStUD[0]->Fill(ASqChNeStUD_east_rp);
	  hASqChPoStLR[0]->Fill(ASqChPoStLR_east_rp);
	  hASqChNeStLR[0]->Fill(ASqChNeStLR_east_rp);
	  hApAmChStUD[0]->Fill(ApAmChStUD_east_rp);
	  hApAmChStLR[0]->Fill(ApAmChStLR_east_rp);
	     
	  ChPlSqUD[0]->Fill(v2ch_east_rp,ASqChPoUD_east_rp);
	  ChPlSqStUD[0]->Fill(v2ch_east_rp,ASqChPoStUD_east_rp);
	     
	  ChMiSqUD[0]->Fill(v2ch_east_rp,ASqChNeUD_east_rp);
	  ChMiSqStUD[0]->Fill(v2ch_east_rp,ASqChNeStUD_east_rp);
	     
	  ChPlMiUD[0]->Fill(v2ch_east_rp,ApAmChUD_east_rp);
	  ChPlMiStUD[0]->Fill(v2ch_east_rp,ApAmChStUD_east_rp);
	     
	  ChPlSqLR[0]->Fill(v2ch_east_rp,ASqChPoLR_east_rp);
	  ChPlSqStLR[0]->Fill(v2ch_east_rp,ASqChPoStLR_east_rp);
	     
	  ChMiSqLR[0]->Fill(v2ch_east_rp,ASqChNeLR_east_rp);
	  ChMiSqStLR[0]->Fill(v2ch_east_rp,ASqChNeStLR_east_rp);
	  
	  ChPlMiLR[0]->Fill(v2ch_east_rp,ApAmChLR_east_rp);
	  ChPlMiStLR[0]->Fill(v2ch_east_rp,ApAmChStLR_east_rp);	    			  
	}

      int iflag2=1;

      if((AChPoUD_west_rp=CalA(NumOfChPoUD_west_rp))>0.8) {iflag2=0;}
      if((AChPoStUD_west_rp=CalA(NumOfChPoStUD_west_rp))>0.8) {iflag2=0;}
      if((AChNeUD_west_rp=CalA(NumOfChNeUD_west_rp))>0.8) {iflag2=0;}
      if((AChNeStUD_west_rp=CalA(NumOfChNeStUD_west_rp))>0.8) {iflag2=0;}
      
      if((AChPoLR_west_rp=CalA(NumOfChPoLR_west_rp))>0.8) {iflag2=0;}
      if((AChPoStLR_west_rp=CalA(NumOfChPoStLR_west_rp))>0.8) {iflag2=0;}
      if((AChNeLR_west_rp=CalA(NumOfChNeLR_west_rp))>0.8) {iflag2=0;}
      if((AChNeStLR_west_rp=CalA(NumOfChNeStLR_west_rp))>0.8) {iflag2=0;}

      if(iflag2==1)
	{
	  ASqChPoUD_west_rp=pow(AChPoUD_west_rp,2);ASqChPoStUD_west_rp=pow(AChPoStUD_west_rp,2);
	  ASqChNeUD_west_rp=pow(AChNeUD_west_rp,2);ASqChNeStUD_west_rp=pow(AChNeStUD_west_rp,2);
	  ApAmChUD_west_rp=AChPoUD_west_rp*AChNeUD_west_rp;ApAmChStUD_west_rp=AChPoStUD_west_rp*AChNeStUD_west_rp;
	  	  
	  ASqChPoLR_west_rp=pow(AChPoLR_west_rp,2);ASqChPoStLR_west_rp=pow(AChPoStLR_west_rp,2);
	  ASqChNeLR_west_rp=pow(AChNeLR_west_rp,2);ASqChNeStLR_west_rp=pow(AChNeStLR_west_rp,2);
	  ApAmChLR_west_rp=AChPoLR_west_rp*AChNeLR_west_rp;ApAmChStLR_west_rp=AChPoStLR_west_rp*AChNeStLR_west_rp;
	     
	  Double_t delta_west_rp=(ASqChPoUD_west_rp-ASqChPoStUD_west_rp+ASqChNeUD_west_rp-ASqChNeStUD_west_rp)/2-(ASqChPoLR_west_rp-ASqChPoStLR_west_rp+ASqChNeLR_west_rp-ASqChNeStLR_west_rp)/2-(ApAmChUD_west_rp-ApAmChStUD_west_rp-ApAmChLR_west_rp+ApAmChStLR_west_rp);

	  hDelta[1]->Fill(delta_west_rp);	     

	  hASqChPoUD[1]->Fill(ASqChPoUD_west_rp);
	  hASqChNeUD[1]->Fill(ASqChNeUD_west_rp);
	  hASqChPoLR[1]->Fill(ASqChPoLR_west_rp);
	  hASqChNeLR[1]->Fill(ASqChNeLR_west_rp);
	  hApAmChUD[1]->Fill(ApAmChUD_west_rp);
	  hApAmChLR[1]->Fill(ApAmChLR_west_rp);
	  
	  hASqChPoStUD[1]->Fill(ASqChPoStUD_west_rp);
	  hASqChNeStUD[1]->Fill(ASqChNeStUD_west_rp);
	  hASqChPoStLR[1]->Fill(ASqChPoStLR_west_rp);
	  hASqChNeStLR[1]->Fill(ASqChNeStLR_west_rp);
	  hApAmChStUD[1]->Fill(ApAmChStUD_west_rp);
	  hApAmChStLR[1]->Fill(ApAmChStLR_west_rp);
	     
	  ChPlSqUD[1]->Fill(v2ch_west_rp,ASqChPoUD_west_rp);
	  ChPlSqStUD[1]->Fill(v2ch_west_rp,ASqChPoStUD_west_rp);
	     
	  ChMiSqUD[1]->Fill(v2ch_west_rp,ASqChNeUD_west_rp);
	  ChMiSqStUD[1]->Fill(v2ch_west_rp,ASqChNeStUD_west_rp);
	     
	  ChPlMiUD[1]->Fill(v2ch_west_rp,ApAmChUD_west_rp);
	  ChPlMiStUD[1]->Fill(v2ch_west_rp,ApAmChStUD_west_rp);
	     
	  ChPlSqLR[1]->Fill(v2ch_west_rp,ASqChPoLR_west_rp);
	  ChPlSqStLR[1]->Fill(v2ch_west_rp,ASqChPoStLR_west_rp);
	     
	  ChMiSqLR[1]->Fill(v2ch_west_rp,ASqChNeLR_west_rp);
	  ChMiSqStLR[1]->Fill(v2ch_west_rp,ASqChNeStLR_west_rp);
	  
	  ChPlMiLR[1]->Fill(v2ch_west_rp,ApAmChLR_west_rp);
	  ChPlMiStLR[1]->Fill(v2ch_west_rp,ApAmChStLR_west_rp);	    			  
	}

      Double_t gamma_pp_east_rp,gamma_nn_east_rp,gamma_pn_east_rp;
      Double_t gamma_pp_west_rp,gamma_nn_west_rp,gamma_pn_west_rp;    
      Double_t gamma_pp_rp,gamma_nn_rp,gamma_pn_rp;

      gamma_pp_east_rp=p_sumcos_east_rp*p_sumcos_east_rp-p_sumcos2_east_rp-p_sumsin_east_rp*p_sumsin_east_rp+p_sumsin2_east_rp;
      gamma_nn_east_rp=n_sumcos_east_rp*n_sumcos_east_rp-n_sumcos2_east_rp-n_sumsin_east_rp*n_sumsin_east_rp+n_sumsin2_east_rp;
      gamma_pn_east_rp=p_sumcos_east_rp*n_sumcos_east_rp-p_sumsin_east_rp*n_sumsin_east_rp;
      
      gamma_pp_east_rp/=p_num_east*(p_num_east-1);
      gamma_nn_east_rp/=n_num_east*(n_num_east-1);
      gamma_pn_east_rp/=p_num_east*n_num_east;

      gamma_pp_west_rp=p_sumcos_west_rp*p_sumcos_west_rp-p_sumcos2_west_rp-p_sumsin_west_rp*p_sumsin_west_rp+p_sumsin2_west_rp;
      gamma_nn_west_rp=n_sumcos_west_rp*n_sumcos_west_rp-n_sumcos2_west_rp-n_sumsin_west_rp*n_sumsin_west_rp+n_sumsin2_west_rp;
      gamma_pn_west_rp=p_sumcos_west_rp*n_sumcos_west_rp-p_sumsin_west_rp*n_sumsin_west_rp;
      
      gamma_pp_west_rp/=p_num_west*(p_num_west-1);
      gamma_nn_west_rp/=n_num_west*(n_num_west-1);
      gamma_pn_west_rp/=p_num_west*n_num_west;
  
      gamma_pp_rp=(p_sumcos_east_rp+p_sumcos_west_rp)*(p_sumcos_east_rp+p_sumcos_west_rp)-(p_sumcos2_east_rp+p_sumcos2_west_rp)-(p_sumsin_east_rp+p_sumsin_west_rp)*(p_sumsin_east_rp+p_sumsin_west_rp)+(p_sumsin2_east_rp+p_sumsin2_west_rp);
      gamma_nn_rp=(n_sumcos_east_rp+n_sumcos_west_rp)*(n_sumcos_east_rp+n_sumcos_west_rp)-(n_sumcos2_east_rp+n_sumcos2_west_rp)-(n_sumsin_east_rp+n_sumsin_west_rp)*(n_sumsin_east_rp+n_sumsin_west_rp)+(n_sumsin2_east_rp+n_sumsin2_west_rp);
      gamma_pn_rp=(p_sumcos_east_rp+p_sumcos_west_rp)*(n_sumcos_east_rp+n_sumcos_west_rp)-(p_sumsin_east_rp+p_sumsin_west_rp)*(n_sumsin_east_rp+n_sumsin_west_rp);

      gamma_pp_rp/=(p_num_east+p_num_west)*(p_num_east+p_num_west-1);
      gamma_nn_rp/=(n_num_east+n_num_west)*(n_num_east+n_num_west-1);
      gamma_pn_rp/=(p_num_east+p_num_west)*(p_num_east+p_num_west);
      

      hpp_east[0]->Fill(gamma_pp_east_rp);
      hpp_west[0]->Fill(gamma_pp_west_rp);
      hnn_east[0]->Fill(gamma_nn_east_rp);
      hnn_west[0]->Fill(gamma_nn_west_rp);
      hpn_east[0]->Fill(gamma_pn_east_rp);
      hpn_west[0]->Fill(gamma_pn_west_rp);
      hpp[0]->Fill(gamma_pp_rp);
      hnn[0]->Fill(gamma_nn_rp);
      hpn[0]->Fill(gamma_pn_rp);    
 
      //**********************************EP************************************

      double AChPoUD_east_ep,AChNeUD_east_ep;
      double ASqChPoUD_east_ep,ASqChNeUD_east_ep;
      double ApAmChUD_east_ep;
      
      double AChPoLR_east_ep,AChNeLR_east_ep;
      double ASqChPoLR_east_ep,ASqChNeLR_east_ep;
      double ApAmChLR_east_ep;
  
      double AChPoUD_west_ep,AChNeUD_west_ep;
      double ASqChPoUD_west_ep,ASqChNeUD_west_ep;
      double ApAmChUD_west_ep;
      
      double AChPoLR_west_ep,AChNeLR_west_ep;
      double ASqChPoLR_west_ep,ASqChNeLR_west_ep;
      double ApAmChLR_west_ep;

      double AChPoStUD_east_ep,AChNeStUD_east_ep;
      double ASqChPoStUD_east_ep,ASqChNeStUD_east_ep;
      double ApAmChStUD_east_ep;
      
      double AChPoStLR_east_ep,AChNeStLR_east_ep;
      double ASqChPoStLR_east_ep,ASqChNeStLR_east_ep;
      double ApAmChStLR_east_ep;
  
      double AChPoStUD_west_ep,AChNeStUD_west_ep;
      double ASqChPoStUD_west_ep,ASqChNeStUD_west_ep;
      double ApAmChStUD_west_ep;
      
      double AChPoStLR_west_ep,AChNeStLR_west_ep;
      double ASqChPoStLR_west_ep,ASqChNeStLR_west_ep;
      double ApAmChStLR_west_ep;   


      Double_t v2ch_east_ep,v2ch_west_ep,v2full_ep;

      v2ch_east_ep=(p_v2sum_east_ep+n_v2sum_east_ep)/(p_num_east+n_num_east);
      v2ch_west_ep=(p_v2sum_west_ep+n_v2sum_west_ep)/(p_num_west+n_num_west);      
      v2full_ep=(p_v2sum_ep+n_v2sum_ep)/(p_num_east+n_num_east+p_num_west+n_num_west);   

      hv2ch[2]->Fill(v2ch_east_ep);
      hv2ch[3]->Fill(v2ch_west_ep);
      hv2[1]->Fill(v2full_ep);     

      //==========================Angle difference between east and west for events with v2==0 ======
      if(fabs(v2ch_east_ep) < 1e-4)
	  hanglediff_eastv2_westep->Fill(EP_east - EP_west); 
      else if(fabs(v2ch_west_ep) < 1e-4)
	  hanglediff_westv2_eastep->Fill(EP_east - EP_west); 
 


      int iflag3=1;
      
      if((AChPoUD_east_ep=CalA(NumOfChPoUD_east_ep))>0.8) {iflag3=0;}
      if((AChPoStUD_east_ep=CalA(NumOfChPoStUD_east_ep))>0.8) {iflag3=0;}
      if((AChNeUD_east_ep=CalA(NumOfChNeUD_east_ep))>0.8) {iflag3=0;}
      if((AChNeStUD_east_ep=CalA(NumOfChNeStUD_east_ep))>0.8) {iflag3=0;}
	 
      if((AChPoLR_east_ep=CalA(NumOfChPoLR_east_ep))>0.8) {iflag3=0;}
      if((AChPoStLR_east_ep=CalA(NumOfChPoStLR_east_ep))>0.8) {iflag3=0;}
      if((AChNeLR_east_ep=CalA(NumOfChNeLR_east_ep))>0.8) {iflag3=0;}
      if((AChNeStLR_east_ep=CalA(NumOfChNeStLR_east_ep))>0.8) {iflag3=0;}


      if(iflag3==1)
	{
	  ASqChPoUD_east_ep=pow(AChPoUD_east_ep,2);ASqChPoStUD_east_ep=pow(AChPoStUD_east_ep,2);
	  ASqChNeUD_east_ep=pow(AChNeUD_east_ep,2);ASqChNeStUD_east_ep=pow(AChNeStUD_east_ep,2);
	  ApAmChUD_east_ep=AChPoUD_east_ep*AChNeUD_east_ep;ApAmChStUD_east_ep=AChPoStUD_east_ep*AChNeStUD_east_ep;
	  	  
	  ASqChPoLR_east_ep=pow(AChPoLR_east_ep,2);ASqChPoStLR_east_ep=pow(AChPoStLR_east_ep,2);
	  ASqChNeLR_east_ep=pow(AChNeLR_east_ep,2);ASqChNeStLR_east_ep=pow(AChNeStLR_east_ep,2);
	  ApAmChLR_east_ep=AChPoLR_east_ep*AChNeLR_east_ep;ApAmChStLR_east_ep=AChPoStLR_east_ep*AChNeStLR_east_ep;
	     
	  Double_t delta_east_ep=(ASqChPoUD_east_ep-ASqChPoStUD_east_ep+ASqChNeUD_east_ep-ASqChNeStUD_east_ep)/2-(ASqChPoLR_east_ep-ASqChPoStLR_east_ep+ASqChNeLR_east_ep-ASqChNeStLR_east_ep)/2-(ApAmChUD_east_ep-ApAmChStUD_east_ep-ApAmChLR_east_ep+ApAmChStLR_east_ep);

	  hDelta[2]->Fill(delta_east_ep);	     

	  hASqChPoUD[2]->Fill(ASqChPoUD_east_ep);
	  hASqChNeUD[2]->Fill(ASqChNeUD_east_ep);
	  hASqChPoLR[2]->Fill(ASqChPoLR_east_ep);
	  hASqChNeLR[2]->Fill(ASqChNeLR_east_ep);
	  hApAmChUD[2]->Fill(ApAmChUD_east_ep);
	  hApAmChLR[2]->Fill(ApAmChLR_east_ep);
	  
	  hASqChPoStUD[2]->Fill(ASqChPoStUD_east_ep);
	  hASqChNeStUD[2]->Fill(ASqChNeStUD_east_ep);
	  hASqChPoStLR[2]->Fill(ASqChPoStLR_east_ep);
	  hASqChNeStLR[2]->Fill(ASqChNeStLR_east_ep);
	  hApAmChStUD[2]->Fill(ApAmChStUD_east_ep);
	  hApAmChStLR[2]->Fill(ApAmChStLR_east_ep);
	     
	  ChPlSqUD[2]->Fill(v2ch_east_ep,ASqChPoUD_east_ep);
	  ChPlSqStUD[2]->Fill(v2ch_east_ep,ASqChPoStUD_east_ep);
	     
	  ChMiSqUD[2]->Fill(v2ch_east_ep,ASqChNeUD_east_ep);
	  ChMiSqStUD[2]->Fill(v2ch_east_ep,ASqChNeStUD_east_ep);
	     
	  ChPlMiUD[2]->Fill(v2ch_east_ep,ApAmChUD_east_ep);
	  ChPlMiStUD[2]->Fill(v2ch_east_ep,ApAmChStUD_east_ep);
	     
	  ChPlSqLR[2]->Fill(v2ch_east_ep,ASqChPoLR_east_ep);
	  ChPlSqStLR[2]->Fill(v2ch_east_ep,ASqChPoStLR_east_ep);
	     
	  ChMiSqLR[2]->Fill(v2ch_east_ep,ASqChNeLR_east_ep);
	  ChMiSqStLR[2]->Fill(v2ch_east_ep,ASqChNeStLR_east_ep);
	  
	  ChPlMiLR[2]->Fill(v2ch_east_ep,ApAmChLR_east_ep);
	  ChPlMiStLR[2]->Fill(v2ch_east_ep,ApAmChStLR_east_ep);	    			  
	}
   
      int iflag4=1;

      if((AChPoUD_west_ep=CalA(NumOfChPoUD_west_ep))>0.8) {iflag4=0;}
      if((AChPoStUD_west_ep=CalA(NumOfChPoStUD_west_ep))>0.8) {iflag4=0;}
      if((AChNeUD_west_ep=CalA(NumOfChNeUD_west_ep))>0.8) {iflag4=0;}
      if((AChNeStUD_west_ep=CalA(NumOfChNeStUD_west_ep))>0.8) {iflag4=0;}
      
      if((AChPoLR_west_ep=CalA(NumOfChPoLR_west_ep))>0.8) {iflag4=0;}
      if((AChPoStLR_west_ep=CalA(NumOfChPoStLR_west_ep))>0.8) {iflag4=0;}
      if((AChNeLR_west_ep=CalA(NumOfChNeLR_west_ep))>0.8) {iflag4=0;}
      if((AChNeStLR_west_ep=CalA(NumOfChNeStLR_west_ep))>0.8) {iflag4=0;}

      if(iflag4==1)
	{
	  ASqChPoUD_west_ep=pow(AChPoUD_west_ep,2);ASqChPoStUD_west_ep=pow(AChPoStUD_west_ep,2);
	  ASqChNeUD_west_ep=pow(AChNeUD_west_ep,2);ASqChNeStUD_west_ep=pow(AChNeStUD_west_ep,2);
	  ApAmChUD_west_ep=AChPoUD_west_ep*AChNeUD_west_ep;ApAmChStUD_west_ep=AChPoStUD_west_ep*AChNeStUD_west_ep;
	  	  
	  ASqChPoLR_west_ep=pow(AChPoLR_west_ep,2);ASqChPoStLR_west_ep=pow(AChPoStLR_west_ep,2);
	  ASqChNeLR_west_ep=pow(AChNeLR_west_ep,2);ASqChNeStLR_west_ep=pow(AChNeStLR_west_ep,2);
	  ApAmChLR_west_ep=AChPoLR_west_ep*AChNeLR_west_ep;ApAmChStLR_west_ep=AChPoStLR_west_ep*AChNeStLR_west_ep;
	     
	  Double_t delta_west_ep=(ASqChPoUD_west_ep-ASqChPoStUD_west_ep+ASqChNeUD_west_ep-ASqChNeStUD_west_ep)/2-(ASqChPoLR_west_ep-ASqChPoStLR_west_ep+ASqChNeLR_west_ep-ASqChNeStLR_west_ep)/2-(ApAmChUD_west_ep-ApAmChStUD_west_ep-ApAmChLR_west_ep+ApAmChStLR_west_ep);

	  hDelta[3]->Fill(delta_west_ep);	     

	  hASqChPoUD[3]->Fill(ASqChPoUD_west_ep);
	  hASqChNeUD[3]->Fill(ASqChNeUD_west_ep);
	  hASqChPoLR[3]->Fill(ASqChPoLR_west_ep);
	  hASqChNeLR[3]->Fill(ASqChNeLR_west_ep);
	  hApAmChUD[3]->Fill(ApAmChUD_west_ep);
	  hApAmChLR[3]->Fill(ApAmChLR_west_ep);
	  
	  hASqChPoStUD[3]->Fill(ASqChPoStUD_west_ep);
	  hASqChNeStUD[3]->Fill(ASqChNeStUD_west_ep);
	  hASqChPoStLR[3]->Fill(ASqChPoStLR_west_ep);
	  hASqChNeStLR[3]->Fill(ASqChNeStLR_west_ep);
	  hApAmChStUD[3]->Fill(ApAmChStUD_west_ep);
	  hApAmChStLR[3]->Fill(ApAmChStLR_west_ep);
	     
	  ChPlSqUD[3]->Fill(v2ch_west_ep,ASqChPoUD_west_ep);
	  ChPlSqStUD[3]->Fill(v2ch_west_ep,ASqChPoStUD_west_ep);
	     
	  ChMiSqUD[3]->Fill(v2ch_west_ep,ASqChNeUD_west_ep);
	  ChMiSqStUD[3]->Fill(v2ch_west_ep,ASqChNeStUD_west_ep);
	     
	  ChPlMiUD[3]->Fill(v2ch_west_ep,ApAmChUD_west_ep);
	  ChPlMiStUD[3]->Fill(v2ch_west_ep,ApAmChStUD_west_ep);
	     
	  ChPlSqLR[3]->Fill(v2ch_west_ep,ASqChPoLR_west_ep);
	  ChPlSqStLR[3]->Fill(v2ch_west_ep,ASqChPoStLR_west_ep);
	     
	  ChMiSqLR[3]->Fill(v2ch_west_ep,ASqChNeLR_west_ep);
	  ChMiSqStLR[3]->Fill(v2ch_west_ep,ASqChNeStLR_west_ep);
	  
	  ChPlMiLR[3]->Fill(v2ch_west_ep,ApAmChLR_west_ep);
	  ChPlMiStLR[3]->Fill(v2ch_west_ep,ApAmChStLR_west_ep);	    			  
	}

      Double_t gamma_pp_east_ep,gamma_nn_east_ep,gamma_pn_east_ep;
      Double_t gamma_pp_west_ep,gamma_nn_west_ep,gamma_pn_west_ep;    
      Double_t gamma_pp_ep,gamma_nn_ep,gamma_pn_ep;

      gamma_pp_east_ep=p_sumcos_east_ep*p_sumcos_east_ep-p_sumcos2_east_ep-p_sumsin_east_ep*p_sumsin_east_ep+p_sumsin2_east_ep;
      gamma_nn_east_ep=n_sumcos_east_ep*n_sumcos_east_ep-n_sumcos2_east_ep-n_sumsin_east_ep*n_sumsin_east_ep+n_sumsin2_east_ep;
      gamma_pn_east_ep=p_sumcos_east_ep*n_sumcos_east_ep-p_sumsin_east_ep*n_sumsin_east_ep;
      
      gamma_pp_east_ep/=p_num_east*(p_num_east-1);
      gamma_nn_east_ep/=n_num_east*(n_num_east-1);
      gamma_pn_east_ep/=p_num_east*n_num_east;

      gamma_pp_west_ep=p_sumcos_west_ep*p_sumcos_west_ep-p_sumcos2_west_ep-p_sumsin_west_ep*p_sumsin_west_ep+p_sumsin2_west_ep;
      gamma_nn_west_ep=n_sumcos_west_ep*n_sumcos_west_ep-n_sumcos2_west_ep-n_sumsin_west_ep*n_sumsin_west_ep+n_sumsin2_west_ep;
      gamma_pn_west_ep=p_sumcos_west_ep*n_sumcos_west_ep-p_sumsin_west_ep*n_sumsin_west_ep;
      
      gamma_pp_west_ep/=p_num_west*(p_num_west-1);
      gamma_nn_west_ep/=n_num_west*(n_num_west-1);
      gamma_pn_west_ep/=p_num_west*n_num_west;
  
      gamma_pp_ep=p_sumcos_ep*p_sumcos_ep-p_sumcos2_ep-p_sumsin_ep*p_sumsin_ep+p_sumsin2_ep;
      gamma_nn_ep=n_sumcos_ep*n_sumcos_ep-n_sumcos2_ep-n_sumsin_ep*n_sumsin_ep+n_sumsin2_ep;
      gamma_pn_ep=p_sumcos_ep*n_sumcos_ep-p_sumsin_ep*n_sumsin_ep;

      gamma_pp_ep/=(p_num_east+p_num_west)*(p_num_east+p_num_west-1);
      gamma_nn_ep/=(n_num_east+n_num_west)*(n_num_east+n_num_west-1);
      gamma_pn_ep/=(p_num_east+p_num_west)*(p_num_east+p_num_west);
         

      hpp_east[1]->Fill(gamma_pp_east_ep);
      hpp_west[1]->Fill(gamma_pp_west_ep);
      hnn_east[1]->Fill(gamma_nn_east_ep);
      hnn_west[1]->Fill(gamma_nn_west_ep);
      hpn_east[1]->Fill(gamma_pn_east_ep);
      hpn_west[1]->Fill(gamma_pn_west_ep);
      hpp[1]->Fill(gamma_pp_ep);
      hnn[1]->Fill(gamma_nn_ep);
      hpn[1]->Fill(gamma_pn_ep);  

      nEvent++; 

      if(nEvent%10000==0)
	cout<<"nEvent = "<<nEvent<<endl;
     
    }
  

  sprintf(tmp,"./root/simulation_%d.root",argc);

  TFile* f = new TFile(tmp,"Recreate");


  hphi_p_east->Write(); 
  hphi_n_east->Write();  
  hphi_p_west->Write();  
  hphi_n_west->Write(); 

  hRes->Write();

  for(int i=0;i<3;i++)
    { 
      hEP[i]->Write();
    }      

  for(int i=0;i<2;i++)
    {
      hpp_east[i]->Write(); 
      hnn_east[i]->Write(); 
      hpn_east[i]->Write(); 
      
      hpp_west[i]->Write();
      hnn_west[i]->Write(); 
      hpn_west[i]->Write(); 
      
      hpp[i]->Write(); 
      hnn[i]->Write(); 
      hpn[i]->Write();

      hv2[i]->Write();
    } 
 
  for(int i=0;i<4;i++)
    {
      
      hv2ch[i]->Write();
      hDelta[i]->Write();
      hASqChPoUD[i]->Write(); 
      hASqChNeUD[i]->Write();
      hASqChPoLR[i]->Write(); 
      hASqChNeLR[i]->Write();
      hApAmChUD[i]->Write();
      hApAmChLR[i]->Write();
      
      hASqChPoStUD[i]->Write();
      hASqChNeStUD[i]->Write();
      hASqChPoStLR[i]->Write();
      hASqChNeStLR[i]->Write();  
      hApAmChStUD[i]->Write();
      hApAmChStLR[i]->Write();
  
      ChPlSqUD[i]->Write();
      ChPlSqStUD[i]->Write();
      
      ChMiSqUD[i]->Write();
      ChMiSqStUD[i]->Write();
      
      ChPlMiUD[i]->Write();
      ChPlMiStUD[i]->Write();
      
      ChPlSqLR[i]->Write();
      ChPlSqStLR[i]->Write();
      
      ChMiSqLR[i]->Write();
      ChMiSqStLR[i]->Write();
      
      ChPlMiLR[i]->Write();
      ChPlMiStLR[i]->Write();
    }

  f->Write();
  f->Close();

}
