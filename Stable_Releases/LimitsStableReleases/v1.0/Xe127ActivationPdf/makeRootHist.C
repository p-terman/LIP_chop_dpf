TH2D* makeRootHist(const char*infile="xe127_ER_PDF_131019d2.tsv", const int nbx=80, const double xmin=0, const double xmax=40, const int nby=250, const double ymin=0.5, const double ymax=3.0){
  TH2D* rootHist = new TH2D("rootHist",Form("%s",infile),nbx,xmin,xmax,nby,ymin,ymax);
  TH2D* hist2To30 = new TH2D("hist2To30",Form("%s",infile),56,2,30,250,0.5,3.0);
  ifstream in;
  in.open(infile);

  float binContent;
  for(int i=1;i<=nbx;i++){
    for(int j=1;j<=nby;j++){
      in>>binContent;
      rootHist->SetBinContent(i,j,binContent);
      if(i>4&&i<61){
	hist2To30->SetBinContent(i-4,j,binContent);
      }
    }
  }
  return hist2To30;
}
