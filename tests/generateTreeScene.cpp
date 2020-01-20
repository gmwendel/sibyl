// Generate Sibyl tests
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom3.h>

#include <RAT/DS/PMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/EV.hh>

#include <string>
#include <iostream>
using namespace std;
using namespace RAT::DS;

PMTInfo* BuildPMTInfo();
void DrawEvent(EV* ev, int time);
void DrawTree(EV* ev, int time);
void DrawCandle(EV* ev, int time);
void DrawSnow(EV* ev, int time);
void hit(EV* ev, int id, double charge);
void FlakeCol(EV* ev, int col, int t0, int time);

TRandom3* rnd = new TRandom3();

int main()
{
  TFile* tfile = new TFile("holiday.root", "RECREATE");
  TTree* runT = new TTree("runT", "runT");

  Run* run = new Run();
  run->SetID(0);
  run->SetType(1);
  run->SetStartTime(1.0);
  run->SetPMTInfo( BuildPMTInfo() );
  runT->Branch("run", run->ClassName(), &run);
  runT->Fill();

  TTree* T = new TTree("T", "T");
  Root* ds = new Root();
  T->Branch("ds", &ds);

  for(int i=0; i<1000; i++)
  {
    ds->PruneEV();
    EV *ev = ds->AddNewEV();
    DrawEvent(ev, i);
    T->Fill();
  }

  tfile->Write();

  return 0;
}

PMTInfo* BuildPMTInfo()
{
  double radius = 6000;
  TVector3 dir(1.0, 0.0, 0.0);
  int type = 1;
  string model = "Santa";

  PMTInfo* pmtinfo = new PMTInfo();
  int nrow = 40;
  int ncol = 100;
  for(int row=0; row<nrow; row++)
  {
    for(int col=0; col<ncol; col++)
    {
      double z = radius - 2*radius/nrow*row;
      double phi = (2*3.1416)/ncol * col;
      double x = radius*cos(phi);
      double y = radius*sin(phi);
      TVector3 pos(x, y, z);
      pmtinfo->AddPMT(pos, dir, type);
    }
  }
  cout << "Count Mid: " << pmtinfo->GetPMTCount();
  int top_start = pmtinfo->GetPMTCount();

  // Top Cap
  int crow = 40;
  for(int row=-1; row<crow+1; row++)
  {
    for(int col=-1; col<crow+1; col++)
    {
      double z = radius+100.0;
      double x = radius - 2*radius*row/crow;
      double y = radius - 2*radius*col/crow;
      TVector3 pos(x, y, z);
      if( sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()) < radius )
        pmtinfo->AddPMT(pos, dir, type);
    }
  }
  int bot_start = pmtinfo->GetPMTCount();
  int top_end = bot_start - 1;
  
  // Bot Cap
  for(int row=-1; row<crow+1; row++)
  {
    for(int col=-1; col<crow+1; col++)
    {
      double z = -(radius+100.0);
      double x = radius - 2*radius*row/crow;
      double y = radius - 2*radius*col/crow;
      TVector3 pos(x, y, z);
      if( sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()) < radius )
        pmtinfo->AddPMT(pos, dir, type);
    }
  }
  int bot_end= pmtinfo->GetPMTCount()-1;

  return pmtinfo;
}

void DrawEvent(EV* ev, int time)
{
  //for(int j=0; j<45; j++)
  //{
  //  PMT* pmt = ev->AddNewPMT();
  //  pmt->SetID(j+time);
  //  pmt->SetTime(1.0);
  //  pmt->SetCharge(j);
  //}
  // Order from back to front
  // Tree
  DrawTree(ev, time);
  DrawCandle(ev, time);
  DrawSnow(ev, time);
  // Snow
}

void DrawTree(EV* ev, int time)
{
  int rows = 40;
  int cols = 100;
  double width = 15.0;
  int middle = 25;
  TRandom3 r;
  
  for(int row=0; row<rows; row++)
  {
    int piece = floor(row*width/rows);
    //hit(ev, row*cols, 1.0);
    for(int c=0; c<cols; c++)
    {
      if( (c > (middle-piece)) && (c < (middle+piece) ) )
      {
        double charge = 0;
        if( !(row%4) && !(c%2) )
          charge = r.Rndm();
        else if( !((row+2)%4) && !((c+1)%2) )
          charge = r.Rndm();
        else
          charge = 0.43;
        hit(ev, row*cols+c, charge); 
      }
    }
    if(row == 1)
    {
      hit(ev, row*cols+middle, 0.700);
    }
    if(row == 2)
    {
      hit(ev, row*cols+middle-2, 0.700);
      hit(ev, row*cols+middle-1, 0.700);
      hit(ev, row*cols+middle, 0.700);
      hit(ev, row*cols+middle+1, 0.700);
      hit(ev, row*cols+middle+2, 0.700);
    }
    if(row == 3)
    {
      hit(ev, row*cols+middle-1, 0.700);
      hit(ev, row*cols+middle, 0.700);
      hit(ev, row*cols+middle+1, 0.700);
    }
    if(row == 4)
    {
      hit(ev, row*cols+middle-2, 0.700);
      hit(ev, row*cols+middle-1, 0.700);
      hit(ev, row*cols+middle, 0.700);
      hit(ev, row*cols+middle+1, 0.700);
      hit(ev, row*cols+middle+2, 0.700);
    }
    if(row == 5)
    {
      hit(ev, row*cols+middle, 0.700);
    }
  }
}

void DrawCandle(EV* ev, int time)
{
  int rows = 40;
  int cols = 100;
  int width = 2.0;
  double gap = 2.0;
  int middle = 75 - width*4 - gap*4;
  double stclr = 0.25;

  for(int stick=0; stick < 9; stick++)
  {
    for(int h=0; h<10; h++) // rows 
    {
      for(int w=0; w<width; w++) // cols 
      {
        hit(ev, (h+22)*cols+ middle + w + stick*(gap+width), stclr);
      }
    }
  } 
  for(int h=0; h<2; h++)
    for(int w=0; w<(width+8*(gap+width)); w++)
      hit(ev, (h+32)*cols + middle + w, stclr);
  for(int h=0; h<(40-22+5); h++)
    for(int w=0; w<width; w++)
      hit(ev, (39-h)*cols + middle+4*(gap+width) + w, stclr);
  // Flames
  double fout = 0.75;
  double fin = 0.82;
  int m=0;
  for(int stick=0; stick < 9; stick++)
  {
    if( stick == 4 )
      m=5;
    else
      m=0;
    // Three layers
    hit(ev, (18-m)*cols + middle + 1 + stick*(gap+width), fout);
    hit(ev, (19-m)*cols + middle + stick*(gap+width), fout);
    hit(ev, (19-m)*cols + middle + 1 + stick*(gap+width), fin);
    hit(ev, (20-m)*cols + middle + 1 + stick*(gap+width), fout);
  }
}

void DrawSnow(EV* ev, int time)
{
  int rows = 40;
  int cols = 100;
  int width = 2.0;
  double gap = 2.0;
  int middle = 75 - width*4 - gap*4;
  double snow= 1.0;
  TRandom3 r;

  for(int k=1; k<10; k++)
  {
    for(int c=0; c<cols; c++)
    {
      int t0 = int(r.Rndm()*200*k);
      FlakeCol(ev, c, t0, time);
    }
  }

  for(int h=0; h<2; h++) // rows 
  {
    for(int w=0; w<100; w++) // cols 
    {
      hit(ev, (39-h)*cols+ w , 1.0);
    }
  }

  // Snow on every pixel, start every 100 frames
}

void FlakeCol(EV* ev, int col, int t0, int time)
{
  int rows = 40;
  int cols = 100;
  if( (time > t0) && ( (time-t0) < 40 ) )
    hit(ev, (cols*(time-t0)) + col, 1.0);
}

void hit(EV* ev, int id, double charge)
{
  PMT* pmt = ev->AddNewPMT();
  pmt->SetID(id);
  pmt->SetTime(1.0);
  pmt->SetCharge(charge);
}
