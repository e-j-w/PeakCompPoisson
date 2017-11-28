FILE *config;
int spectrum[NSPECT], startCh[NSPECT], endCh[NSPECT], numSpectra, endSpectrum,
    maxNumCh, numSimData;
double sfc[NSPECT], sfw[NSPECT];
int addBackground; // 0=no, 1=constant background
int plotOutput;    // 0=no, 1=yes, 2=detailed
int saveStats;     // 0=no, 1=yes
int saveBG;        // 0=no, 1=yes
int saveResults;   // 0=no, 1=yes
char expDataName[256], simDataName[256], statsDataName[256],
    resultsDataName[256], bgDataName[256]; // file names
char str[256], str1[256], str2[256];

void readConfigFile(const char *fileName) {
  int index = 0;
  numSpectra = 0;
  endSpectrum = 0;
  maxNumCh = 0;
  numSimData = 0;
  if ((config = fopen(fileName, "r")) == NULL) {
    printf("ERROR: Cannot open the config file %s!\n", fileName);
    exit(-1);
  }
  while (!(feof(config))) // go until the end of file is reached
  {
    if (fgets(str, 256, config) != NULL) {
      if (index < NSPECT)
        // spectrum, channel range and step function parameter data
        if (sscanf(str, "%i %i %i %lf %lf", &spectrum[index], &startCh[index],
                   &endCh[index], &sfc[index], &sfw[index]) == 5) {
          if (spectrum[index] > endSpectrum)
            endSpectrum = spectrum[index];
          if ((endCh[index] - startCh[index] + 1) > maxNumCh)
            maxNumCh = endCh[index] - startCh[index] + 1;
          index++;
          numSpectra++;
        }
      if (sscanf(str, "%s %s", str1, str2) == 2) // single parameter data
      {
        if (strcmp(str1, "EXPERIMENT_DATA") == 0)
          strcpy(expDataName, str2);
        if (strcmp(str1, "SIMULATED_DATA") == 0){
          strcpy(simDataName, str2);
          numSimData++;
        }  
        if (strcmp(str1, "ADD_BACKGROUND") == 0) {
          if (strcmp(str2, "yes") == 0)
            addBackground = 1;
          else
            addBackground = 0;
        }
        if (strcmp(str1, "PLOT_OUTPUT") == 0) {
          if (strcmp(str2, "yes") == 0)
            plotOutput = 1;
          else if (strcmp(str2, "detailed") == 0)
            plotOutput = 2;
          else
            plotOutput = 0;
        }
        if (strcmp(str1, "VERBOSE") == 0) {
          if (strcmp(str2, "yes") == 0)
            verbosity = 1;
          else if (strcmp(str2, "debug") == 0)
            verbosity = 2;
          else
            verbosity = 0;
        }
        if (strcmp(str1, "SAVE_STATS") == 0) {
          if (strcmp(str2, "yes") == 0)
            saveStats = 1;
          else
            saveStats = 0;
        }
        if (strcmp(str1, "STATS_DATA_NAME") == 0)
          strcpy(statsDataName, str2);

        if (strcmp(str1, "SAVE_RESULTS") == 0) {
          if (strcmp(str2, "yes") == 0)
            saveResults = 1;
          else
            saveResults = 0;
        }
        if (strcmp(str1, "SAVE_BACKGROUND") == 0) {
          if (strcmp(str2, "yes") == 0)
            saveBG = 1;
          else
            saveBG = 0;
        }
        if (strcmp(str1, "BACKGROUND_DATA_NAME") == 0)
          strcpy(bgDataName, str2);
      }
      if (sscanf(str, "%s %s", str1, str2) == 1) // listing of simulated data
      {
        if (strcmp(str1, "<---END_OF_PARAMETERS--->") == 0)
          break;
        /*else if (strcmp(str1, "SIMULATED_DATA") != 0) {
          strcpy(simDataName, str1);
          numSimData++;
        }*/
      }
    }
  }
  fclose(config);

  if(verbosity>0){
    printf("Taking experiment data from file: %s\n", expDataName);
    if (numSimData > 0)
      printf("Taking simulated data from file: %s\n", simDataName);
    else {
      printf("ERROR: no simulated data specified!\n");
      exit(-1);
    }
    for (index = 0; index < numSpectra; index++)
      printf("Will compare spectrum %i from channels %i to %i.\n",
            spectrum[index], startCh[index], endCh[index]);
    for (index = 0; index < numSpectra; index++)
      printf("Will fit step function for spectrum %i using centroid %.2lf and "
            "width %.2lf\n",
            spectrum[index], sfc[index], sfw[index]);
    if (addBackground == 0)
      printf("Will not add background to simulated data.\n");
    if (addBackground == 1)
      printf("Will add a linear background to simulated data.\n");
    if (plotOutput == 0)
      printf("Will not plot output data.\n");
    if (plotOutput == 1)
      printf("Will plot output data.\n");
    if (plotOutput == 2)
      printf("Will plot detailed output data.\n");
    if (saveStats == 1)
      printf("Will save fit stats to file %s\n", statsDataName);
    if (saveResults == 1)
      printf("Saving scaled fit results.\n");
    if (saveBG == 1)
      printf("Will save fit background to file %s\n", bgDataName);
    printf("Finished reading parameter file...\n");
  }
  
}
