Double_t resp1(Double_t *x, Double_t *par) {
  return peak1g->Eval(x[0])*par[0];
}

Double_t resp2(Double_t *x, Double_t *par) {
  return peak2g->Eval(x[0])*par[0];
}

Double_t resp3(Double_t *x, Double_t *par) {
  return peak3g->Eval(x[0])*par[0];
}

Double_t expf(Double_t *x, Double_t *par) {
  return TMath::Exp( par[0]+par[1]*x[0] );
}

//Double_t ex_respf(Double_t *x, Double_t *par) {
//  return resp1(x,par)+resp2(x,par+1)+resp2(x,par+2)+expf(x,par+3);
//}
Double_t ex_respf(Double_t *x, Double_t *par) {
  return resp1(x,par)+expf(x,par+1);
}


//---------------------
//For the second peak:
Double_t expf2(Double_t *x, Double_t *par) {
  return TMath::Exp( par[0]+par[1]*x[0] );
}

Double_t resp4(Double_t *x, Double_t *par) {
  return peak4g->Eval(x[0])*par[0];
}

Double_t ex_respf2(Double_t *x, Double_t *par) {
  return resp4(x,par)+expf2(x,par+1);
}
