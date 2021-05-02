#! /usr/local/bin/root.exe
{
   gROOT->LoadMacro("./ShogunReconstructor.C+");
   ShogunReconstructor();
   gROOT->ProcessLine(".q");
}
