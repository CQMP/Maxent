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

#pragma once
#include <alps/hdf5/archive.hpp>
#include <alps/utilities/cast.hpp>
#include <Eigen/Core>
#include <algorithm>

//TODO: this really should be in ALPS, not here.
namespace alps {
    namespace hdf5 {
        template<> struct is_continuous<Eigen::VectorXd>
            : public is_continuous<Eigen::VectorXd::Scalar>
        {};
        template<> struct is_continuous<Eigen::VectorXd const >
            : public is_continuous<Eigen::VectorXd::Scalar>
        {};
        namespace detail {
          template<> struct get_extent<Eigen::VectorXd> {
                  static std::vector<std::size_t> apply(Eigen::VectorXd const & value) {
                      using alps::hdf5::get_extent;
                      std::vector<std::size_t> result(1, value.size());
                      if (value.size()) {
                          std::vector<std::size_t> first(get_extent(value[0]));
                          if (true)
                              for(std::size_t i=1;i<value.size();i++) {
                                  std::vector<std::size_t> size(get_extent(value[i]));
                                  if (
                                         first.size() != size.size()
                                      || !std::equal(first.begin(), first.end(), size.begin())
                                  )
                                      throw archive_error("no rectangular matrix" + ALPS_STACKTRACE);
                              }
                          std::copy(first.begin(), first.end(), std::back_inserter(result));
                      }
                      return result;
                  }
          };
        }
        void save(
            archive &ar
          , std::string const & path
          , Eigen::VectorXd &value
          , std::vector<std::size_t> size = std::vector<std::size_t>()
          , std::vector<std::size_t> chunk =std::vector<std::size_t>()
          , std::vector<std::size_t> offset = std::vector<std::size_t>()
       );
    }
}
