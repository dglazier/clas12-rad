#pragma once

namespace rad{
  namespace clas12{
    constexpr short UndefPDG=1E4;

    //regions
    constexpr short FT = 1000;
    constexpr short FD = 2000;
    constexpr short CD = 3000;
    constexpr short BD = 4000;

  
    //detectors
    constexpr short FTOF = 12;
    constexpr short CTOF = 4;
    constexpr short CND  = 3;
    constexpr short CVT   = 5;
    constexpr short DC   = 6;
    constexpr short EC   = 7;
    constexpr short ECAL   = 7;
    constexpr short FTCAL   = 10;
    constexpr short FTTRK   = 13;
    constexpr short FTHODO   = 11;
    constexpr short HTCC   = 15;
    constexpr short LTCC   = 16;
    constexpr short BMT   = 1;
    constexpr short FMT   = 8;
    constexpr short RF   = 17;
    constexpr short RICH   = 18;
    constexpr short RTPC   = 19;
    constexpr short HEL   = 20;
    constexpr short BAND   = 21;

    class DetId2Name {

      std::map<short,string> _map;

    public :
      DetId2Name(){
	_map[FTOF]="FTOF";
	_map[CTOF]="CTOF";
	_map[CND]="CND";
	_map[CVT]="CVT";
	_map[DC]="DC";
	_map[ECAL]="ECAL";
	_map[FTCAL]="FTCAL";
	_map[FTTRK]="FTTRK";
	_map[FTHODO]="FTHODO";
	_map[HTCC]="HTCC";
	_map[LTCC]="LTCC";
	_map[BMT]="BMT";
	_map[FMT]="FMT";
	_map[RF]="RF";
	_map[RICH]="RICH";
	_map[RTPC]="RTPC";
	_map[HEL]="HEL";
	_map[BAND]="BAND";
      }
      const string& DetName(short idet){return _map[idet];}
    };
    // _map[]="";

    //layers
    constexpr short FTOF1A = 1;
    constexpr short FTOF1B = 2;
    constexpr short FTOF2 = 3;
    //CDET scint layers same as detectors
    //constexpr short CND  = 3;
    //constexpr short CTOF = 4;
 
    constexpr short CNDOFF = 150; //CND1-CNDOFF == actual layer number, stil stops conflict with FTOF scintillator
    constexpr short CND1  = 151;
    constexpr short CND2  = 152;
    constexpr short CND3  = 153;
 
    constexpr short PCAL   = 1;
    constexpr short ECIN   = 4;
    constexpr short ECOUT  = 7;

    //additional TRAJECTORY layers
    constexpr short DC1  = 6;
    constexpr short DC2  = 12;
    constexpr short DC3  = 18;
    constexpr short DC4  = 24;
    constexpr short DC5  = 30;
    constexpr short DC6  = 36;

    constexpr short CVT1  = 1;
    constexpr short CVT2  = 2;
    constexpr short CVT3  = 3;
    constexpr short CVT4  = 4;
    constexpr short CVT5  = 5;
    constexpr short CVT6  = 6;
    constexpr short CVT7  = 7;
    constexpr short CVT8  = 8;
    constexpr short CVT9  = 9;
    constexpr short CVT10  = 10;
    constexpr short CVT11  = 11;
    constexpr short CVT12  = 12;

    constexpr short BANDOFF = 250;
    constexpr short BTOF1 = 251;
    constexpr short BTOF2 = 252;
    constexpr short BTOF3 = 253;
    constexpr short BTOF4 = 254;
    constexpr short BTOF5 = 255;
    constexpr short BVETO = 256;

 
  }
}
