/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#include <iostream>
#pragma once

void write_minimal_param_file(const std::string &str){
    std::ofstream tmpfile(str.c_str());
    tmpfile<<"BETA=2" <<std::endl;
    tmpfile<<"X_0=.1" <<std::endl;
    tmpfile<<"X_1=.2" <<std::endl;
    tmpfile<<"X_2=.3" <<std::endl;
    tmpfile<<"X_3=.4" <<std::endl;
    tmpfile<<"SIGMA_0=.5" <<std::endl;
    tmpfile<<"SIGMA_1=.5" <<std::endl;
    tmpfile<<"SIGMA_2=.5" <<std::endl;
    tmpfile<<"SIGMA_3=.5" <<std::endl;
    tmpfile.close();
}
void write_minimal_input_file(const std::string &str){
  std::ofstream tmpfile(str.c_str());
  tmpfile<<0<<" "<<0.1<<" "<<0.5<<std::endl;
  tmpfile<<1<<" "<<0.2<<" "<<0.5<<std::endl;
  tmpfile<<2<<" "<<0.3<<" "<<0.5<<std::endl;
  tmpfile<<3<<" "<<0.4<<" "<<0.5<<std::endl;
  tmpfile<<4<<" "<<0.5<<" "<<0.5<<std::endl;
  tmpfile.close();

}
//default model 
void write_minimal_tab_file(const std::string &str){
  {
    std::ofstream tmpfile(str.c_str());
    tmpfile<<-20<<" "<<1./40<<std::endl;
    tmpfile<<20<<" "<<1./40<<std::endl;
    tmpfile.close();
  }
}
void write_tab_file_with_junk(const std::string &str){
  {
    std::ofstream tmpfile(str.c_str());
    tmpfile<<"# this is a file with junk in it."<<std::endl;
    tmpfile<<-20<<" "<<1./40<<" default model data" <<std::endl;
    tmpfile<<20<<" "<<1./40<<std::endl;
    tmpfile.close();
  }
}
