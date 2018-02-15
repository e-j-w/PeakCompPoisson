#include "peak_comp_poisson_root.h"
#include "read_config.c"

int main(int argc, char *argv[]) {
  // Creating an instance of TApplication
  // This is evidently needed for auto-loading of ROOT libraries,
  // otherwise the program may crash upon execution depending on how ROOT
  // is set up.
  // http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=14064
  int ac;
  char *av[10];
  theApp = new TApplication("App", &ac, av);

  int i = 0;
  int j = 0;
  int k = 0;
  FILE *expData, *simData, *results, *scalingFactors;

  if (argc != 2) {
    printf("\npeak_comp_poisson parameter_file\n");
    printf("Compares the .mca spectra designated in the parameter file "
           "specified and generates cool statistics.\n\n");
    exit(-1);
  }

  // initialize values
  addBackground = 0;
  for (i = 0; i < NSPECT; i++)
    for (j = 0; j < S32K; j++) {
      expHist[i][j] = 0.;
      for (k = 0; k < NSIMDATA; k++)
        {
          simHist[k][i][j] = 0.;
          scaledSimHist[k][i][j] = 0.;
        }
        
    }

  for (i = 0; i < NPAR; i++)
    for (j = 0; j < NSPECT; j++)
      aFinal[i][j] = 0.;

  chisq = 0.;

  // initialize ROOT stuff
  char simName[132];
  char dataName[132];
  for (i = 0; i < NSPECT; i++) {
    sprintf(simName, "sim_%2d", i);
    sim[i] = new TH1D(simName, ";;", S32K, 0, S32K - 1);
  }
  for (i = 0; i < NSPECT; i++) {
    sprintf(dataName, "data_%2d", i);
    data[i] = new TH1D(dataName, ";;", S32K, 0, S32K - 1);
  }

  readConfigFile(argv[1]); // grab data from the config file

  // check that the number of spectra being compared is fine
  if (endSpectrum >= NSPECT) {
    printf("ERROR: A spectrum number specified in the parameter file is larger "
           "than the maximum value of %i.  Reduce it or increase NSPECT in "
           "peak_comp.h and recompile.\n",
           NSPECT);
    exit(-1);
  }

  // read in the .mca files
  // check file extension of exp data and copy into float histogram
  if ((expData = fopen(expDataName, "r")) == NULL) {
    printf("ERROR: Cannot open the experiment data file %s!\n", expDataName);
    exit(-1);
  }
  const char *dot = strrchr(expDataName, '.'); // get the file extension
  if (strcmp(dot + 1, "mca") == 0)
    readMCA(expData, expDataName, expHist);
  else if (strcmp(dot + 1, "fmca") == 0)
    readFMCA(expData, expDataName, expHist);
  else {
    printf("ERROR: Improper type of input file: %s\n", expDataName);
    printf(
        "Integer array (.mca) and float array (.fmca) files are supported.\n");
    exit(-1);
  }
  fclose(expData);

  for (i = 1; i < NSPECT; i++)
    for (j = 0; j < S32K; j++)
      data[i]->Fill(j, expHist[i][j]);
  
  for(i = 0; i < numSimData; i++)
    {
      if ((simData = fopen(simDataName[i], "r")) == NULL) {
        printf("ERROR: Cannot open the simulated data file %s!\n", simDataName[i]);
        exit(-1);
      }
      const char *dots = strrchr(simDataName[i], '.'); // get the file extension
      if (strcmp(dots + 1, "mca") == 0)
        readMCA(simData, simDataName[i], simHist[i]);
      else if (strcmp(dots + 1, "fmca") == 0)
        readFMCA(simData, simDataName[i], simHist[i]);
      else {
        printf("ERROR: Improper type of input file: %s\n", simDataName[i]);
        printf(
            "Integer array (.mca) and float array (.fmca) files are supported.\n");
        exit(-1);
      }
      fclose(simData);
    }

  /*// read into ROOT
  for (i = 0; i < numSpectra; i++)
    for (j = 0; j < S32K; j++)
      sim[i]->Fill(j, simHist[0][i][j]);*/

  if(verbosity>0)
    printf("Spectra read in...\n");

  
  find_chisqMin();

  // scale simulated data
  for (i = 0; i < numSpectra; i++)
    for (j = 0; j < S32K; j++)
      for(k=0;k<numSimData;k++)
        {
          scaledSimHist[k][spectrum[i]][j]=0.;
          scaledSimHist[k][spectrum[i]][j] += aFinal[k+3][i] * simHist[k][spectrum[i]][j];
        }

  // add background to simulated data
  for (i = 0; i < numSpectra; i++)
    for (j = 0; j < S32K; j++){
      if (addBackground == 1){
        bgHist[spectrum[i]][j] =
          aFinal[0][i] +
          aFinal[1][i] * erfc(((double)j - sfc[i]) / sfw[i]);
      }else if (addBackground == 2){
        bgHist[spectrum[i]][j] =
          aFinal[0][i] +
          aFinal[1][i] * (double)j;
      }else if (addBackground == 3){
        bgHist[spectrum[i]][j] = aFinal[0][i];
      }else{
        bgHist[spectrum[i]][j] = 0.;
      }
      
    }
      

  // fit result histogram
  for (i = 0; i < numSpectra; i++)
    for (j = 0; j < S32K; j++)
      {
        resultsHist[spectrum[i]][j] = bgHist[spectrum[i]][j];
        for(k=0;k<numSimData;k++)
          resultsHist[spectrum[i]][j] += (float)scaledSimHist[k][spectrum[i]][j];
      }
      

  // print output
  if(verbosity>0)
    printf("Fit chisq: %.15f\n", chisq);
  else
    printf("%.15f\n", chisq);

  // write results
  sprintf(str, "fit_poisson.fmca");
  if ((results = fopen(str, "w")) == NULL) {
    printf("ERROR: Cannot open the output file %s!\n", str);
    exit(-1);
  }
  for (i = 0; i < numSpectra; i++) {
    fwrite(resultsHist[spectrum[i]], S32K * sizeof(float), 1, results);
  }
  fclose(results);

  // save scaling factors
  sprintf(str, "scalingFactors.dat");
  if ((scalingFactors = fopen(str, "w")) == NULL) {
    printf("ERROR: Cannot open the output file %s!\n", str);
    exit(-1);
  }
  for (i = 0; i < numSpectra; i++) {
    fprintf(scalingFactors, "%d %.9f\n", i + 1, aFinal[0][i]);
  }
  fclose(scalingFactors);

  /*printf("Scaling factors:\n");
  for (i = 0; i < numSpectra; i++) {
    printf("%d %.9f\n", i + 1, aFinal[0][i]);
  }*/

  // plot results
  if (plotOutput >= 1)
    plotSpectra();

  return 0; // great success
}

double lrchisq(const double *par) {
  // chisq from likelihood ratio test
  // for more information see Baker and Cousins pg. 439 and Appendix
  double yi = 0.;      // model
  double ni = 0.;      // experiment
  double lrchisq = 0.; // likelihood ratio chisq
  int i = 0;
  int j = 0;

  for (i = startCh[spCurrent]; i <= endCh[spCurrent]; i++) {
    ni = expCurrent[i]; // events in ith bin

    yi=0;
    for(j=0; j<numSimData; j++)
      {
        // calculate model in the ith bin
        yi += par[j+3] * simCurrent[j][i];
      }
    // add background if neccesary
    if (addBackground == 1){
      yi += par[0] + par[1] * erfc(((double)i - sfc[spCurrent]) / sfw[spCurrent]);
    }else if (addBackground == 2){
      yi += par[0] + par[1] * (double)i;
    }else if (addBackground == 3){
      yi += par[0];
    }

    // evaluate chisq given input parameters
    if (ni > 0.)
      lrchisq += (yi - ni + ni * log(ni / yi));
    else
      lrchisq += yi; // the log(0) case
  }
  lrchisq *= 2.;

  /*printf("Parameters: %f %f %f %f, computed lrchisq: %f\n",par[0],par[1],par[2],par[3],lrchisq);
  getc(stdin);*/
  return lrchisq;
}

double pchisq(const double *par) {
  // pearson chisq
  // for more information see Baker and Cousins pg. 438 and Appendix
  double yi = 0.;     // model
  double ni = 0.;     // experiment
  double pchisq = 0.; // pearson chisq
  int i = 0;
  int j = 0;

  for (i = startCh[spCurrent]; i <= endCh[spCurrent]; i++) {
    ni = expCurrent[i]; // events in ith bin

    if (ni < 0.)
      ni = 0.;

    yi=0;
    for(j=0; j<numSimData; j++)
      {
        // calculate model in the ith bin
        yi += par[j+3] * simCurrent[j][i];
      }
    // add background if neccesary
    if (addBackground == 1){
      yi += par[0] + par[1] * erfc(((double)i - sfc[spCurrent]) / sfw[spCurrent]);
    }else if (addBackground == 2){
      yi += par[0] + par[1] * (double)i;
    }else if (addBackground == 3){
      yi += par[0];
    }
    

    // evaluate chisq given input parameters
    pchisq += (ni - yi) * (ni - yi) / yi;
  }
  return pchisq;
}

double nchisq(const double *par) {
  // neyman chisq
  // for more information see Baker and Cousins pg. 438 and Appendix
  double yi = 0.;     // model
  double ni = 0.;     // experiment
  double nchisq = 0.; // neyman ratio chisq
  int i = 0;
  int j = 0;

  for (i = startCh[spCurrent]; i <= endCh[spCurrent]; i++) {
    ni = expCurrent[i]; // events in ith bin

    yi=0;
    for(j=0; j<numSimData; j++)
      {
        // calculate model in the ith bin
        yi += par[j+3] * simCurrent[j][i];
      }
    // add background if neccesary
    if (addBackground == 1){
      yi += par[0] + par[1] * erfc(((double)i - sfc[spCurrent]) / sfw[spCurrent]);
    }else if (addBackground == 2){
      yi += par[0] + par[1] * (double)i;
    }else if (addBackground == 3){
      yi += par[0];
    }

    // evaluate chisq given input parameters
    if (ni > 0.) // have to drop all bins where ni = 0
      nchisq += (ni - yi) * (ni - yi) / ni;
  }
  return nchisq;
}

void find_chisqMin() {
  int i = 0;
  int j = 0;
  int k = 0;
  char str[256];
  double intExp,intSim;

  if(verbosity>0)
    printf("Fitting data...\n");

  for (i = 0; i < numSpectra; i++) {
    // for more information see minimizer class documentation
    // https://root.cern.ch/root/html/ROOT__Math__Minimizer.html
    char minName[132] = "Minuit";
    char algoName[132] = ""; // default (Migard for Minuit)
    ROOT::Math::Minimizer *min =
        ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    // set tolerance , etc...
    min->SetMaxFunctionCalls(10000000); // for Minuit
    min->SetMaxIterations(100000);
    min->SetTolerance(0.001);
    if(verbosity>1)
      min->SetPrintLevel(1);
    else
      min->SetPrintLevel(0); // set to 1 for more info

    // communication to the likelihood ratio function
    // (via global parameters)
    for (j = 0; j < S32K; j++) {
      expCurrent[j] = (double)expHist[spectrum[i]][j];
    }
    for (j = 0; j < S32K; j++) {
      for (k=0;k<numSimData;k++)
        simCurrent[k][j] = (double)simHist[k][spectrum[i]][j];
    }
    spCurrent = i;

    //calculate integrals
    intSim=0.;
    intExp=0.;
    for (j = startCh[spCurrent]; j <= endCh[spCurrent]; j++){
      for (k=0;k<numSimData;k++)
        intSim += (double)simHist[k][spectrum[i]][j];
      intExp += (double)expHist[spectrum[i]][j]; 
    }

    /*//print simulated values fed into minimizer
    double expSum=0.,simSum=0.;
    for (j = startCh[spCurrent]; j <= endCh[spCurrent]; j++){
      printf("expCurrent[%i] = %f   simCurrent[%i] = %f\n",j,expCurrent[j],j,simCurrent[j]);
      expSum+=expCurrent[j];
      simSum+=simCurrent[j];  
    }
    printf("sum exp = %f   sum sim = %f, exp/sim ratio = %f\n",expSum,simSum,expSum/simSum);
    getc(stdin);*/

    // create function wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor lr(&lrchisq, 3+NSIMDATA); // likelihood ratio chisq
    /* ROOT::Math::Functor lr(&nchisq,3+NSIMDATA); // neyman chisq */
    /* ROOT::Math::Functor lr(&pchisq,3+NSIMDATA); // pearson chisq - inf problems! */

    // step size and starting variables
    // may need to change for best performance
    // under different running conditions
    double ratio = intExp/intSim;
    double variable[3+NSIMDATA];
    double step[3+NSIMDATA];

    for (j=0;j<3+NSIMDATA;j++){
      variable[j] = ratio/2.;
      step[j] = ratio;
    }

    // 94Sr (low stats)
    // double step[3] = {0.00001,0.0001,0.0001};
    // double variable[3] = {0.0001,0.001,0.001};

    // 84Kr (high stats)
    // double step[3] = {0.001,0.001,0.001};
    // double variable[3] = {0.01,0.01,0.01};

    min->SetFunction(lr);

    // Set pars for minimization
    
    for (j=0;j<3+NSIMDATA;j++){
      sprintf(str, "a%i", j);
      min->SetVariable(j, str, variable[j], step[j]);
    }
    
    /*min->SetVariable(0, "a0", variable[0], step[0]);
    min->SetVariable(1, "a1", variable[1], step[1]);
    min->SetVariable(2, "a2", variable[2], step[2]);*/

    // variable limits (optional)
    /* min->SetVariableLimits(0,0.,1E3); */
    /* min->SetVariableLimits(1,0.,1E3); */
    /* min->SetVariableLimits(2,0.,1E3); */

    // do the minimization
    min->Minimize();

    // grab parameters from minimum
    const double *xs = min->X();

    // print results
    // degrees of freedom assuming 3 pars
    /*double ndf = (endCh[spCurrent] - startCh[spCurrent] + 1) - 3;
    const double *exs = min->Errors();
    std::cout << "spectrum " << i << " " << xs[0] << " " << exs[0] << " "
              << xs[1] << " " << exs[1] << " " << xs[2] << " " << exs[2] << " "
              << "chi2: " << min->MinValue()
              << " chi2/ndf: " << min->MinValue() / ndf << std::endl;
    std::cout << i << " chi2: " << min->MinValue()
              << " chi2/ndf: " << min->MinValue() / ndf << std::endl;
    std::cout << min->MinValue() << " ";*/

    if(verbosity>0)
      printf("Spectrum %i fit parameters: %f %f %f\n",i,xs[0],xs[1],xs[2]);

    // assuming 3 parameters
    // save pars
    for (j = 0; j < 3+NSIMDATA; j++)
      aFinal[j][i] = xs[j];

    // add to total chisq
    /* if(i==0 || i==4) // groups 1 and 5 only */
    /* if(i==1 || i==3) // groups 2 and 4 only */
    chisq += min->MinValue();
  }
}

void plotSpectra() {
  int i, j, k;
  TH1D *results[NSPECT];
  TH1D *resultsBGData[NSPECT];
  TH1D *resultsSimData[NSIMDATA][NSPECT];
  TCanvas *c = new TCanvas("c1", "", 1618, 1000);

  // initialize and fill results histo
  char resultsName[132],resultsBGName[132],resultsSimDataName[132];
  for (i = 0; i < NSPECT; i++) {
    sprintf(resultsName, "results_%2d", i);
    results[i] = new TH1D(resultsName, ";;", S32K, 0, S32K - 1);
    if(plotMode>=1)
      {
        sprintf(resultsBGName, "resultsBG_%2d", i);
        resultsBGData[i] = new TH1D(resultsName, ";;", S32K, 0, S32K - 1);
        for(k=0;k<numSimData;k++)
          {
            sprintf(resultsSimDataName, "resultsSimData_%2d_%2d", i, k);
            resultsSimData[k][i] = new TH1D(resultsSimDataName, ";;", S32K, 0, S32K - 1);
          }
      }
    
  }
  // be careful with indicies here
  for (i = 0; i < numSpectra; i++)
    for (j = 0; j < S32K; j++)
      {
        results[spectrum[i]]->Fill(j, resultsHist[spectrum[i]][j]);
        if(plotMode>=1)
          {
            resultsBGData[spectrum[i]]->Fill(j, bgHist[spectrum[i]][j]);
            for(k=0;k<numSimData;k++)
              resultsSimData[k][spectrum[i]]->Fill(j, scaledSimHist[k][spectrum[i]][j]);
          }
      }

  // display limits
  Double_t low[NSPECT], high[NSPECT];

  // divide canvas to accomodate all spectra
  // assumes number to plot <= 3*2=6
  c->Divide(3, 2, 1E-6, 1E-6, 0);

  // formatting the plot area
  // be careful with indicies
  for (i = 1; i <= 6; i++) {
    c->GetPad(i)->SetRightMargin(0.01);
    c->GetPad(i)->SetTopMargin(0.025);
    Int_t ind = i - 1;
    low[i] = startCh[ind];
    high[i] = endCh[ind];
  }

  // plot
  for (i = 1; i <= 6; i++) {
    c->cd(i);
    // data in black
    data[i]->SetLineStyle(1);
    data[i]->SetLineWidth(2);
    data[i]->SetLineColor(12);
    data[i]->GetXaxis()->SetRangeUser(low[i], high[i]);
    data[i]->SetStats(0);
    data[i]->Draw();

    // simulation in red
    results[i]->SetLineStyle(1);
    results[i]->SetLineWidth(2);
    results[i]->SetLineColor(46);
    results[i]->Draw("SAME");
    if(plotMode>=1)
      {
        //plot individual simulated data
        for(k=0;k<numSimData;k++)
          {
            resultsSimData[k][i]->SetLineStyle(1);
            resultsSimData[k][i]->SetLineWidth(1);
            resultsSimData[k][i]->SetLineColor(799+k*20);
            resultsSimData[k][i]->Draw("SAME");
          }
        if(plotMode!=2)
          {
            //plot background
            resultsBGData[i]->SetLineStyle(1);
            resultsBGData[i]->SetLineWidth(1);
            resultsBGData[i]->SetLineColor(920);
            resultsBGData[i]->Draw("SAME");
          }
        
      }
    
  }

  c->cd(1);
  // x1,y1,x2,y2
  TLegend *leg = new TLegend(0.70, 0.85, 0.90, 0.95);
  leg->AddEntry(data[1], "Experiment", "l");
  leg->AddEntry(results[1], "Simulation", "l");
  if(plotMode>=1)
    {
      for(k=0;k<numSimData;k++)
        {
          sprintf(resultsName, "Branch %i", k+1);
          leg->AddEntry(resultsSimData[k][1], resultsName, "l");
        }
      
    }
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->Draw();

  c->SetBorderMode(0);
  c->Update();

  theApp->Run(kTRUE);
}

int readMCA(FILE *inp, char *filename, float inpHist[NSPECT][S32K]) {
  int i = 0;
  int j = 0;

  int mcaHist[NSPECT][S32K];
  for (i = 0; i <= endSpectrum; i++)
    if (fread(mcaHist[i], S32K * sizeof(int), 1, inp) != 1) {
      printf("ERROR: Error reading file %s!\n", filename);
      exit(-1);
    }

  for (i = 0; i < NSPECT; i++)
    for (j = 0; j < S32K; j++)
      inpHist[i][j] = (float)mcaHist[i][j];

  return 1;
}

int readFMCA(FILE *inp, char *filename, float inpHist[NSPECT][S32K]) {
  int i = 0;

  for (i = 0; i <= endSpectrum; i++)
    if (fread(inpHist[i], S32K * sizeof(float), 1, inp) != 1) {
      printf("ERROR: Error reading file %s!\n", filename);
      exit(-1);
    }

  return 1;
}
