{
    gROOT->LoadMacro("simulation_v2.C++");
    int i = 550;
    while(i > 0){
        simulation_v2(1208, i, 200000);
        i -= 100;
    }
}
