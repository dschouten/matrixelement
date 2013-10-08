//
// author Doug Schouten <doug dot schouten at triumf dot ca>
//

#include "matrix/Fortran.hh"
#include "matrix/Constants.hh"

namespace hepstd 
{
  
  void prepareCommonBlocks( double MH, double MT, double WH )
  {
    static bool printedCON  = false;
    static bool printedSM   = false;
    static bool printedHEFT = false;
    static bool printedMASS = false;
    static bool printedRS   = false;
    
    const Complex COMPLEXI	= Complex(0,1)			;
    
    const double AEWM1		= 1.3250698e+02			;
    const double GF		= 1.1663900e-05			; 
    const double CKM22		= 1.000   			;
    const double AS		= 0.118				;
    const double YMB		= 4.200				;
    const double YMT		= 164.5				;
    const double YMTAU		= 1.777				;
    const double LRS            = 3.000e+03                     ;
    
    const double SQRT__AS	= sqrt(AS)			;
    const double G		= 2.0*SQRT__AS*sqrt(M_PI)	;
    const double MZ__EXP__2	= pow(zMass,2)			;
    const double MZ__EXP__4	= pow(zMass,4)			;
    const double AEW		= 1.0/AEWM1			;
    const double SQRT__AEW	= sqrt(AEW)			;
    const double EE		= 2.0*SQRT__AEW*sqrt(M_PI)	;
    const double MW__EXP__2	= pow(wMass,2)			;
    const double SW2		= 1.0-MW__EXP__2/MZ__EXP__2	;
    const double CW		= sqrt(1.0-SW2)			;
    const double SQRT__SW2	= sqrt(SW2)			;
    const double SW		= SQRT__SW2			;
    const double G1		= EE/CW				;
    const double GW		= EE/SW				;
    const double V		= (2.0*wMass*SW)/EE		;
    const double EE__EXP__2	= pow(EE,2)			;
    const double MW__EXP__12	= pow(wMass,12)			;
    const double MW__EXP__10	= pow(wMass,10)			;
    const double MW__EXP__8	= pow(wMass,8)			;
    const double MW__EXP__6	= pow(wMass,6)			;
    const double MW__EXP__4	= pow(wMass,4)			;
    const double V__EXP__2	= pow(V,2)			;
    const double YB		= (YMB*M_SQRT_2)/V		;
    const double YT		= (YMT*M_SQRT_2)/V		;
    const double YTAU		= (YMTAU*M_SQRT_2)/V		;
    const double GW__EXP__2	= pow(GW,2)			;
    const double CW__EXP__2	= pow(CW,2)			;
    const double SW__EXP__2	= pow(SW,2)			;
    const double G__EXP__2	= pow(G,2)			;
    const double KAPPA          = 2.0/LRS                       ;
    
    const double MH__EXP__2		= pow(MH,2)			;
    const double MH__EXP__4		= pow(MH,4)			;
    const double MH__EXP__6		= pow(MH,6)			;
    const double MH__EXP__8		= pow(MH,8)			;
    const double MH__EXP__10		= pow(MH,10)			;
    const double MH__EXP__12		= pow(MH,12)			;
    const double MT__EXP__2		= pow(MT,2)			;
    const double MT__EXP__4		= pow(MT,4)			;
    const double MT__EXP__6		= pow(MT,6)			;   
    const double LAM			= MH__EXP__2/(2.0*V__EXP__2)	;
    const double MUH			= sqrt(LAM*V__EXP__2)		;
    const long double R__EXP__6              = pow( MH / MT, 6 )             ;
    const long double R__EXP__4              = pow( MH / MT, 4 )             ;
    const long double R__EXP__2              = pow( MH / MT, 2 )             ;
    const long double GH		     = ( -(G__EXP__2*(1.000000E+00+(1.300000E+01*R__EXP__6)/(1.680000E+04)+R__EXP__4/(1.680000E+02)
							      +(7.000000E+00*R__EXP__2)/(1.200000E+02)))/(1.200000E+01*M_PI_2*V) ) / M_2PI;
    
    // const double GPHI			= ( -(G__EXP__2*(1.000000E+00+MH__EXP__6/
    // 						 (5.600000E+02
    // 						  *MT__EXP__6)+MH__EXP__4/(9.000000E+01*MT__EXP__4)+MH__EXP__2/
    // 						 (1.200000E+01*MT__EXP__2)))/(8.000000E+00*M_PI_2*V) );
    
    if( GlobalFlags::debug && !printedCON )
    {
      printedCON = true;
      std::cout << std::endl 
		<< "\t----------------------------- Model Constants ----------------------------- " << std::endl    
		<< "\tCOMPLEXI = " << COMPLEXI << std::endl
		<< "\tAEWM1 = " << AEWM1 << std::endl	
		<< "\tGF = " << GF << std::endl	
		<< "\tCKM22 = " << CKM22 << std::endl	
		<< "\tAS = " << AS << std::endl	
		<< "\tYMB = " << YMB << std::endl	
		<< "\tYMT = " << YMT << std::endl	
		<< "\tYMTAU = " << YMTAU << std::endl	
		<< "\tSQRT__AS = " << SQRT__AS << std::endl
		<< "\tG = " << G << std::endl	
		<< "\tMZ__EXP__2 = " << MZ__EXP__2 << std::endl
		<< "\tMZ__EXP__4 = " << MZ__EXP__4 << std::endl
		<< "\tAEW = " << AEW << std::endl	
		<< "\tSQRT__AEW = " << SQRT__AEW << std::endl
		<< "\tEE = " << EE << std::endl	
		<< "\tMW__EXP__2 = " << MW__EXP__2 << std::endl
		<< "\tSW2 = " << SW2 << std::endl	
		<< "\tCW = " << CW << std::endl	
		<< "\tSQRT__SW2 = " << SQRT__SW2 << std::endl
		<< "\tSW = " << SW << std::endl	
		<< "\tG1 = " << G1 << std::endl	
		<< "\tGW = " << GW << std::endl	
		<< "\tV = " << V << std::endl	
		<< "\tEE__EXP__2 = " << EE__EXP__2 << std::endl
		<< "\tMW__EXP__12 = " << MW__EXP__12 << std::endl
		<< "\tMW__EXP__10 = " << MW__EXP__10 << std::endl
		<< "\tMW__EXP__8 = " << MW__EXP__8 << std::endl
		<< "\tMW__EXP__6 = " << MW__EXP__6 << std::endl
		<< "\tMW__EXP__4 = " << MW__EXP__4 << std::endl
		<< "\tV__EXP__2 = " << V__EXP__2 << std::endl
		<< "\tYB = " << YB << std::endl	
		<< "\tYT = " << YT << std::endl	
		<< "\tYTAU = " << YTAU << std::endl	
		<< "\tGW__EXP__2 = " << GW__EXP__2 << std::endl
		<< "\tCW__EXP__2 = " << CW__EXP__2 << std::endl
		<< "\tSW__EXP__2 = " << SW__EXP__2 << std::endl
		<< "\tG__EXP__2 = " << G__EXP__2 << std::endl
		<< "\tMH__EXP__2 = " << MH__EXP__2 << std::endl	
		<< "\tMH__EXP__4 = " << MH__EXP__4 << std::endl	
		<< "\tMH__EXP__6 = " << MH__EXP__6 << std::endl	
		<< "\tMH__EXP__8 = " << MH__EXP__8 << std::endl	
		<< "\tMH__EXP__10 = " << MH__EXP__10 << std::endl	
		<< "\tMH__EXP__12 = " << MH__EXP__12 << std::endl	
		<< "\tMT__EXP__2 = " << MT__EXP__2 << std::endl	
		<< "\tMT__EXP__4 = " << MT__EXP__4 << std::endl	
		<< "\tMT__EXP__6 = " << MT__EXP__6 << std::endl	
		<< "\tLAM = " << LAM << std::endl		
		<< "\tMUH = " << MUH << std::endl		
		<< "\tR__EXP__6 = " << R__EXP__6 << std::endl              
		<< "\tR__EXP__4 = " << R__EXP__4 << std::endl              
		<< "\tR__EXP__2 = " << R__EXP__2 << std::endl              
		<< "\tGH = " << GH << std::endl;
    }
    
    C_to_mirror.imirror = 1;
    
    C_polarization.pol[0] = 1;
    C_polarization.pol[1] = 1;

    C_tomecomb.iselectedhel = -1;
    
    C_to_matrix.isum_hel = 1024;
    C_to_matrix.multi_channel = 0;
    
    C_masses.MH = MH;
    C_masses.MT = MT;
    C_masses.MB = hepstd::bMass;
    C_masses.MW = hepstd::wMass;
    C_masses.MZ = hepstd::zMass;
    
    C_widths.WH = WH;
    C_widths.WT = hepstd::topDecayWidth( hepstd::tMass );  
    C_widths.WW = hepstd::wWidth;
    C_widths.WZ = hepstd::zWidth;
    
    C_smcouplings.GC_3=Complex(-(EE*COMPLEXI))();
    C_smcouplings.GC_16=Complex((CKM22*EE*COMPLEXI)/(SW*M_SQRT_2))();
    C_smcouplings.GC_1=Complex(-(EE*COMPLEXI)/3.000000e+00)();
    C_smcouplings.GC_7=Complex(CW*COMPLEXI*GW)();
    C_smcouplings.GC_5=Complex(COMPLEXI*G)();
    C_smcouplings.GC_4=Complex(-G)();
    C_smcouplings.GC_31=Complex((EE__EXP__2*COMPLEXI*V)/(2.000000*SW__EXP__2))();
    C_smcouplings.GC_33=Complex(-((COMPLEXI*YB)/M_SQRT_2))();
    C_smcouplings.GC_2=Complex((2.000000e+00*EE*COMPLEXI)/3.000000e+00)();
    C_smcouplings.GC_22=Complex((CW*EE*COMPLEXI)/(2.000000e+00*SW))();
    C_smcouplings.GC_23=Complex(-(EE*COMPLEXI*SW)/(6.000000e+00*CW))();
    C_smcouplings.GC_21=Complex(-(CW*EE*COMPLEXI)/(2.000000e+00*SW))();
    C_smcouplings.GC_24=Complex((EE*COMPLEXI*SW)/(2.000000e+00*CW))();
    C_smcouplings.GC_25=Complex(COMPLEXI*GW*SW)();    
    
    C_heftcouplings.GC_17=Complex((EE*COMPLEXI)/(SW*M_SQRT_2))();
    C_heftcouplings.GC_6=Complex(COMPLEXI*G)();
    C_heftcouplings.GC_5=Complex(-G)();
    C_heftcouplings.GC_9=Complex(-(G*(double)GH))();
    C_heftcouplings.GC_8=Complex(-(COMPLEXI*( (double)GH )))();
    C_heftcouplings.GC_37=Complex((EE__EXP__2*COMPLEXI*V)/(2.000000e+00*SW__EXP__2))();
    
    C_rscouplings.GC_26=Complex((EE*COMPLEXI)/(SW*M_SQRT_2))();
    C_rscouplings.GC_68=Complex(-(EE__EXP__2*COMPLEXI*KAPPA*V__EXP__2)/(8.000000))();
    C_rscouplings.GC_11=Complex(-(COMPLEXI*KAPPA)/2.000000e+00)();
    
    if( GlobalFlags::debug && !printedSM )
    {    	
      std::cout << std::endl
		<< "\t----------------------------- COMMON/SMCOUPLINGS/ ---------------------------- " << std::endl
		<< "\tGC_1 : " << "(" << C_smcouplings.GC_1.real  << "," << C_smcouplings.GC_1.imag  << ")" << std::endl
		<< "\tGC_3 : " << "(" << C_smcouplings.GC_3.real  << "," << C_smcouplings.GC_3.imag  << ")" << std::endl
		<< "\tGC_16: " << "(" << C_smcouplings.GC_16.real << "," << C_smcouplings.GC_16.imag << ")" << std::endl 
		<< "\tGC_7 : " << "(" << C_smcouplings.GC_7.real  << "," << C_smcouplings.GC_7.imag  << ")" << std::endl
		<< "\tGC_5 : " << "(" << C_smcouplings.GC_5.real  << "," << C_smcouplings.GC_5.imag  << ")" << std::endl
		<< "\tGC_4 : " << "(" << C_smcouplings.GC_4.real  << "," << C_smcouplings.GC_4.imag  << ")" << std::endl 
		<< "\tGC_31: " << "(" << C_smcouplings.GC_31.real << "," << C_smcouplings.GC_31.imag << ")" << std::endl
		<< "\tGC_33: " << "(" << C_smcouplings.GC_33.real << "," << C_smcouplings.GC_33.imag << ")" << std::endl
		<< "\tGC_2 : " << "(" << C_smcouplings.GC_2.real  << "," << C_smcouplings.GC_2.imag  << ")" << std::endl 
		<< "\tGC_21: " << "(" << C_smcouplings.GC_21.real << "," << C_smcouplings.GC_21.imag << ")" << std::endl
		<< "\tGC_22: " << "(" << C_smcouplings.GC_22.real << "," << C_smcouplings.GC_22.imag << ")" << std::endl
		<< "\tGC_23: " << "(" << C_smcouplings.GC_23.real << "," << C_smcouplings.GC_23.imag << ")" << std::endl
		<< "\tGC_25: " << "(" << C_smcouplings.GC_25.real << "," << C_smcouplings.GC_25.imag << ")" << std::endl;
      printedSM = true;
    }
    
    if( GlobalFlags::debug && !printedHEFT )
    {
      std::cout << std::endl
		<< "\t----------------------------- COMMON/HEFTCOUPLINGS/ ----------------------------- " << std::endl
		<< "\tGC_5 : " << "(" << C_heftcouplings.GC_5.real  << "," << C_heftcouplings.GC_5.imag  << ")" << std::endl
		<< "\tGC_6 : " << "(" << C_heftcouplings.GC_6.real  << "," << C_heftcouplings.GC_6.imag  << ")" << std::endl
		<< "\tGC_8 : " << "(" << C_heftcouplings.GC_8.real  << "," << C_heftcouplings.GC_8.imag  << ")" << std::endl 
		<< "\tGC_9 : " << "(" << C_heftcouplings.GC_9.real  << "," << C_heftcouplings.GC_9.imag  << ")" << std::endl
		<< "\tGC_17: " << "(" << C_heftcouplings.GC_17.real << "," << C_heftcouplings.GC_17.imag << ")" << std::endl 
		<< "\tGC_37: " << "(" << C_heftcouplings.GC_37.real << "," << C_heftcouplings.GC_37.imag << ")" << std::endl;    
      printedHEFT = true;
    }
    
    if( GlobalFlags::debug && !printedRS )
    {
      std::cout << std::endl
		<< "\t----------------------------- COMMON/RSCOUPLINGS/ ----------------------------- " << std::endl
		<< "\tGC_11: " << "(" << C_rscouplings.GC_11.real << "," << C_rscouplings.GC_11.imag << ")" << std::endl
		<< "\tGC_26: " << "(" << C_rscouplings.GC_26.real << "," << C_rscouplings.GC_26.imag << ")" << std::endl    
		<< "\tGC_68: " << "(" << C_rscouplings.GC_68.real << "," << C_rscouplings.GC_68.imag << ")" << std::endl;    
      printedRS = true;
    }
    
    if( GlobalFlags::debug && !printedMASS )
    {
      std::cout << std::endl
		<< "\t---------------------------------- COMMON/MASSES/ --------------------------------- " << std::endl
		<< "\tMB: " << C_masses.MB << std::endl
		<< "\tMH: " << C_masses.MH << std::endl  
		<< "\tMT: " << C_masses.MT << std::endl
		<< "\tMW: " << C_masses.MW << std::endl
		<< "\tMZ: " << C_masses.MZ << std::endl;
      std::cout << std::endl
		<< "\t---------------------------------- COMMON/WIDTHS/ --------------------------------- " << std::endl
		<< "\tWH: " << C_widths.WH << std::endl  
		<< "\tWT: " << C_widths.WT << std::endl
		<< "\tWW: " << C_widths.WW << std::endl
		<< "\tWZ: " << C_widths.WZ << std::endl;
      printedMASS = true;
    }
    
  }
}
