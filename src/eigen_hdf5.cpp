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

#include"eigen_hdf5.hpp"

//TODO: this really should be somewhere in ALPS. Also: why copy size, chunk, and offset?
namespace alps {
    namespace hdf5 {
        void save(
            archive &ar
          , std::string const & path
          , Eigen::VectorXd &value
          , std::vector<std::size_t> size
          , std::vector<std::size_t> chunk
          , std::vector<std::size_t> offset
       ) {
            using alps::cast;
            if (ar.is_group(path))
                ar.delete_group(path);
            if (is_continuous<Eigen::VectorXd>::value && value.size() == 0)
                ar.write(path, static_cast<Eigen::VectorXd::Scalar const *>(NULL), std::vector<std::size_t>());
            else if (is_continuous<Eigen::VectorXd>::value) {
                std::vector<std::size_t> extent(get_extent(value));
                std::copy(extent.begin(), extent.end(), std::back_inserter(size));
                std::copy(extent.begin(), extent.end(), std::back_inserter(chunk));
                std::fill_n(std::back_inserter(offset), extent.size(), 0);
                ar.write(path, value.data(), size, chunk, offset);

            } else if (value.size() == 0)
                ar.write(path, static_cast<int const *>(NULL), std::vector<std::size_t>());
        }
    }
}
