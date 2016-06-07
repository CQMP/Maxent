/*
 * Copyright (C) 1998-2016 ALPS Collaboration
 * 
 *     This program is free software; you can redistribute it and/or modify it
 *     under the terms of the GNU General Public License as published by the Free
 *     Software Foundation; either version 2 of the License, or (at your option)
 *     any later version.
 * 
 *     This program is distributed in the hope that it will be useful, but WITHOUT
 *     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
 *     more details.
 * 
 *     You should have received a copy of the GNU General Public License along
 *     with this program; if not, write to the Free Software Foundation, Inc., 59
 *     Temple Place, Suite 330, Boston, MA 02111-1307 USA
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
