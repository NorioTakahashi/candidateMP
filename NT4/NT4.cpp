#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <NT4.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
ad_comm::change_datafile_name("NT4_ctrl.dat"); // set parameter values for MP
  ssb_pop_target.allocate("ssb_pop_target");
  t_pop.allocate("t_pop");
  tune_year.allocate("tune_year");
  k1_cpue.allocate("k1_cpue");
  k2_cpue.allocate("k2_cpue");
  t1_cpue.allocate("t1_cpue");
  k3_cpue.allocate("k3_cpue");
  k4_cpue.allocate("k4_cpue");
  t2_cpue.allocate("t2_cpue");
  k1_gt_limit.allocate("k1_gt_limit");
  t_gt_limit.allocate("t_gt_limit");
  n_age2_limit.allocate("n_age2_limit");
  swit_cpue.allocate("swit_cpue");
  swit_gt_limit.allocate("swit_gt_limit");
  max_change_up.allocate("max_change_up");
  max_change_down.allocate("max_change_down");
  min_change.allocate("min_change");
  debug_write.allocate("debug_write");
first_cpue_yr = 1969;
last_cpue_yr = 2016;
  hist_cpue.allocate(first_cpue_yr,last_cpue_yr,"hist_cpue");
ad_comm::change_datafile_name("sbtOMdata");
  first_yr.allocate("first_yr");
  last_yr.allocate("last_yr");
  implementation_yr.allocate("implementation_yr");
  current_yr.allocate("current_yr");
  quota.allocate(first_yr-1,implementation_yr-1,"quota");
  proj_cpuen.allocate(first_yr,current_yr-2,"proj_cpuen");
  proj_cpueb.allocate(first_yr,current_yr-2,"proj_cpueb");
  proj_asurv.allocate(first_yr,current_yr-1,"proj_asurv");
  proj_gtn.allocate(first_yr,last_yr,"proj_gtn");
  proj_cvgtn.allocate(first_yr,last_yr,"proj_cvgtn");
  proj_gtrec.allocate(first_yr,last_yr,"proj_gtrec");
ad_comm::change_datafile_name("POP_nsim.dat");
  nPOPdat.allocate("nPOPdat");
ad_comm::change_datafile_name("POP_proj.dat");
  POPproj.allocate(1,nPOPdat,1,5,"POPproj");
ad_comm::change_datafile_name("POP_hist.dat");
  nPOPdat2.allocate("nPOPdat2");
  POPhist.allocate(1,nPOPdat2,1,5,"POPhist");
  obs_cpue.allocate(first_cpue_yr,current_yr-2);
  obs_gtn.allocate(2016,current_yr-2);
  obs_cvgtn.allocate(2016,current_yr-2);
  obs_gtrec.allocate(2016,current_yr-2);
  yrs.allocate(first_cpue_yr,current_yr-2);
  POPind.allocate(2002,current_yr-5);
  Mja.allocate(2002,last_yr-5,2006,last_yr);
  Rja.allocate(2002,last_yr-5,2006,last_yr);
  wja.allocate(2002,last_yr-5,2006,last_yr);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  // primary variables
  double  mu_pop;     // average POP index over most recent years
  double  t_cpue;      // work: number of years over which slope of CPUEs(age4+) is estimated
  double  slope_cpue;  // work: slope of cpue for age4+
  double  k_cpue;      // work: gain parameter for cpue slope
  double  mu_gt;       // average GT estimate over most recent years
  
  int  flg_gt_limit;  // flag of whether GT limit is used
  
  double muCK; // f/ Rich's code
  
  double  TAC_cpue;     // TAC by cpue slope
  double  TAC_gt_limit; // TAC by GT limit
  double  tmp_tac;      // temporary TAC before constraints
  double  tac_change;   // amount of tac change
  double  TAC;         // final specified TAC for current year
  
  // variables (f/ Rich's code)
  int a,c,cc,i,y;
  double wsum;
  // (f/ Rich's code)
  POPind.initialize();
  Mja.initialize();
  Rja.initialize();
  wja.initialize();
  // concatenate the CPUE and GT index data
  for (int ii = first_cpue_yr; ii <= current_yr - 2; ii++) {
  
    yrs(ii) = ii;
    
    if(ii <= last_cpue_yr) {
      obs_cpue(ii) = hist_cpue(ii);  //age4+
    }
    
    if(ii > last_cpue_yr) {
      obs_cpue(ii) = proj_cpuen(ii);  //age4+
    }
    
    // copied Rich's example MP code
    if(ii == 2016){
      obs_gtn(ii) = 2417786.; // hard-wired (norio corrected according to Ann's e-mail 180531)
      obs_cvgtn(ii) = 0.2132; // hard-wired (norio added f/ estimates Excel file)
      obs_gtrec(ii) = 22.;     // hard-wired
    }
    
    if(ii > 2016){
      obs_gtn(ii) = proj_gtn(ii);
      obs_cvgtn(ii) = proj_cvgtn(ii);
      obs_gtrec(ii) = proj_gtrec(ii);
    }
  }
  // ******** CPUE(age4+) slope part of HCR ********
  // CK POP-index (f/ Rich's code)
  
  int cmin = 2002;
  int cmax = current_yr-5;
  int ymin = 2006;
  int ymax = current_yr;
  
  ///////////////
  // POP index //
  ///////////////
  
  // aggregate across adult capture ages
  for(int n=1;n<=nPOPdat2;n++) {
     c = POPhist(n,1);
     y = POPhist(n,2);
     a = POPhist(n,3);
     Mja(c,y) += POPhist(n,5);
     Rja(c,y) += POPhist(n,4);
  }
  for(int n=1;n<=nPOPdat;n++) {
     c = POPproj(n,1);
     y = POPproj(n,2);
     a = POPproj(n,3);
     Mja(c,y) += POPproj(n,5);
     Rja(c,y) += POPproj(n,4);
  }
  // create the weighting for the index and sum the index
  for(c=cmin;c<=cmax;c++) {
     wsum = sum(Rja(c));
     for(y=ymin;y<=ymax;y++) wja(c,y) = Rja(c,y)/wsum;
     for(y=ymin;y<=ymax;y++) {
       if(Rja(c,y) > 0) POPind(c) += (Mja(c,y)/Rja(c,y)) * wja(c,y);
       
     }
     
  }
  //cout << POPind << endl;
  //Iref = sum(POPind(c1,c2))/double(c2-c1+1);
  muCK = 0.;
  //for(y=cmax-tau3+1;y<=cmax;y++) muCK += POPind(y)/double(tau3);
  for(y=cmax-t_pop+1 ;y<=cmax;y++) muCK += POPind(y)/double(t_pop);
  //cout << current_yr << " " << Iref << " " << muCK << " " << delCK << " " << delCPUE << " " << delGT << endl;
  //cout << current_yr << " " << muCK << " " << endl;
  mu_pop = muCK;
  // ---- calculate average CPUE over most recent years
  
  int yr1, yr2;
  if( mu_pop <= ssb_pop_target ){
  
    yr1 = current_yr-2 - t1_cpue + 1;
    yr2 = current_yr-2;
  
    slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /
                 sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );
    k_cpue = slope_cpue < 0 ? k1_cpue : k2_cpue;
  
  }
  else if( mu_pop > ssb_pop_target ){
  
    if(current_yr > tune_year){
      yr1 = current_yr-2 - t2_cpue + 1;
      yr2 = current_yr-2;
    }
    else{
      yr1 = current_yr-2 - t1_cpue + 1;
      yr2 = current_yr-2;
    }
  
    slope_cpue = sum( elem_prod( log(obs_cpue(yr1,yr2))-mean(log(obs_cpue(yr1,yr2))), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) ) /
                 sum( elem_prod( yrs(yr1,yr2)-mean(yrs(yr1,yr2)), yrs(yr1,yr2)-mean(yrs(yr1,yr2)) ) );
    if(current_yr > tune_year){
      k_cpue = slope_cpue < 0 ? k3_cpue : k4_cpue;
    }
    else{
      k_cpue = slope_cpue < 0 ? k1_cpue : k2_cpue;
    }
  
  }
  TAC_cpue = quota(implementation_yr-1) * (1. + k_cpue * slope_cpue);
  // ******** GT(age2) part of HCR ********
  
  // (1) used as limit
  
  //----calculate ratio of (recent mean GT estimate)/(limit value of N_age2)
  // first, need to remove data of when no esitimate is available ("-11" is embedded)
  // the following code was based on the code (C1GT2.tpl) provided by Rich
  yr1 = current_yr - t_gt_limit - 1 < first_yr-2 ? first_yr-2 : current_yr - t_gt_limit - 1;
  yr2 = current_yr-2;
  yr1 = yr1 < 2016 ? 2016 : yr1;
  int cnt;
  cnt = 0;
  
  mu_gt = 0.;
  
  for(int ii=yr1; ii<=yr2; ii++) {
    mu_gt = obs_gtn(ii) != -11. ? mu_gt+obs_gtn(ii) : mu_gt;
    cnt = obs_gtn(ii) != -11. ? cnt+1 : cnt;
  }
  mu_gt /= cnt;
  
  TAC_gt_limit = quota(implementation_yr-1) * k1_gt_limit * (mu_gt/n_age2_limit)*(mu_gt/n_age2_limit);
  
  flg_gt_limit = mu_gt < n_age2_limit ? 1 : 0;
  // (2) used as slope
  // ---- calculate GT slope
  // added based on the code (C1GT2.tpl) provided by Rich
  // ----set TAC
  // use TAC from CPUE slope part of HCR only
  if(swit_cpue==1 && swit_gt_limit==0) tmp_tac = TAC_cpue;
  // use TAC from GT limit part of HCR only
  if(swit_gt_limit==1 && swit_cpue==0){
  
    tmp_tac = flg_gt_limit==1 ? TAC_gt_limit : quota(implementation_yr-1);
  
  }
  // use TAC from CPUE slope and GT limit parts of HCR
  if(swit_cpue==1 && swit_gt_limit==1){
  
    tmp_tac = flg_gt_limit==1 ? min(TAC_gt_limit, TAC_cpue) : TAC_cpue;
    
  }
  
  
  tac_change = tmp_tac - quota(implementation_yr-1);
  
  if( fabs(tac_change) <= min_change ){
  
    TAC = quota(implementation_yr-1);
  
  }
  else if( tac_change > max_change_up ) {
  
    TAC = quota(implementation_yr-1) + max_change_up;
  
  }
  else if( tac_change < -1 * max_change_down) {
  
    TAC = quota(implementation_yr-1) - max_change_down;
  
  }
  else{
  
    TAC = tmp_tac;
  
  }
  
  
  //----output for sbtproj.exe (copied from Rich's test code
  ofstream ofs("tacfile");
  ofs << implementation_yr << " " << TAC << endl;
  //----debug write?
  if(debug_write){
    DebugWrite1(slope_cpue, mu_gt, mu_pop);
  }
  exit(1);
}

void model_parameters::userfunction(void)
{
  f =0.0;
  f=0.;
}

void model_parameters::DebugWrite1(double slope_cpue, double mu_gt, double mu_pop)
{
  ofstream fdebug;
  int ii;
  fdebug.open("debug_write_ctrl.txt", ios::trunc);
  fdebug << "---- check NT*_ctrl.dat inputs ----" << endl;
  fdebug << "ssb_pop_target = " << ssb_pop_target << endl;
  fdebug << "t_pop          = " << t_pop          << endl;
  fdebug << "tune_year      = " << tune_year      << endl;
  fdebug << "k1_cpue        = " << k1_cpue        << endl;
  fdebug << "k2_cpue        = " << k2_cpue        << endl;
  fdebug << "t1_cpue        = " << t1_cpue        << endl;
  fdebug << "k3_cpue        = " << k3_cpue        << endl;
  fdebug << "k4_cpue        = " << k4_cpue        << endl;
  fdebug << "t2_cpue        = " << t2_cpue        << endl;
  fdebug << "k1_gt_limit    = " << k1_gt_limit    << endl;
  fdebug << "t_gt_limit     = " << t_gt_limit     << endl;
  fdebug << "n_age2_limit   = " << n_age2_limit   << endl;
  fdebug << "swit_cpue     = " << swit_cpue     << endl;
  fdebug << "swit_gt_limit = " << swit_gt_limit << endl;
  fdebug << "max_change_up   = " << max_change_up   << endl;
  fdebug << "max_change_down = " << max_change_down << endl;
  fdebug << "min_change      = " << min_change      << endl;
  fdebug << "debug_write  = " << debug_write  << endl;
  for(ii = first_cpue_yr; ii <= last_cpue_yr; ii++){
    fdebug << "hist_cpue(" << ii << ") = " << hist_cpue(ii) << endl;
  }
  fdebug.close();
  fdebug.open("debug_write_cpue_gt.txt", ios::trunc);
  fdebug << "---- check values of cpue, GT, GT CV, POP set to variables ----" << endl;
  for(ii = first_cpue_yr; ii <= current_yr-2; ii++){
    fdebug << "obs_cpue(" << ii << ") = " << obs_cpue(ii) << endl;
  }
  fdebug << endl;
  for(ii = 2016; ii <= current_yr-2; ii++){
    fdebug << "obs_gtn(" << ii << ") = " << obs_gtn(ii) << endl;
  }
  fdebug << endl;
  for(ii = 2016; ii <= current_yr-2; ii++){
    fdebug << "obs_cvgtn(" << ii << ") = " << obs_cvgtn(ii) << endl;
  }
  fdebug << endl;
  for(ii = 2002; ii <= current_yr-5; ii++){
    fdebug << "POPind(" << ii << ") = " << POPind(ii) << endl;
  }
  fdebug.close();
  fdebug.open("debug_write_primary1.txt", ios::app);
  fdebug << "---- check values set to primary variables - 1 ----" << endl;
  fdebug << "current_yr = " << current_yr << endl;
  fdebug << "  slope_cpue = " << slope_cpue << endl;
  fdebug << "  mu_gt      = " << mu_gt      << endl;
  fdebug << "  mu_pop     = " << mu_pop     << endl;
  fdebug << endl;
  fdebug.close();
  return;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(const dvector& gradients){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
