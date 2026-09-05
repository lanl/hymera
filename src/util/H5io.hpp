#pragma once
/// \file   H5io.hpp
/// \brief  Functions for reading and writing scalar and vector data to HDF5
/// files
/// \author Oleksii Beznosov, LANL
/// \date   2024-07-26
/// \note   This file is part of the kinetic code base and is subject to the
/// license terms in the LICENSE file found in the top-level directory of this
/// distribution.

#include <H5Cpp.h>
#include <cassert> // For std::ceil
#include <cmath>   // For std::ceil
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>
#include <Kokkos_Core.hpp>
#include <string>
#include <type_traits>

/// Expand the H5 namespace to include functions for reading and writing scalar
/// and vector data
namespace H5 {
const std::string scalarsGroupName = "/scalars";
const std::string vectorsGroupName = "/vectors";
const std::string viewsGroupName = "/views";

/// Write a scalar value to an HDF5 group
/// \tparam T The type of the scalar value
/// \param group The HDF5 group to write to
/// \param name The name of the attribute
/// \param value The scalar value to write
template <typename T>
void writeScalar(H5::Group &group, const std::string &name, const T &value) {
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
  if constexpr (std::is_same_v<T, int>) {
    datatype = H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, double>) {
    datatype = H5::PredType::NATIVE_DOUBLE;
  } // Add more types as needed

  H5::DataSpace dataspace(H5S_SCALAR);
  H5::Attribute attribute = group.createAttribute(name, datatype, dataspace);
  attribute.write(datatype, &value);
}

/// Write a scalar value to an HDF5 file
/// \tparam T The type of the scalar value
/// \param filename The HDF5 file to write to, defoulting to group
/// scalarGroupName
/// \param name The name of the attribute
/// \param value The scalar value to write
template <typename T>
void writeScalar(const std::string &filename, const std::string &name,
                 const T &value) {
  H5::H5File file(filename, H5F_ACC_RDWR);
  H5::Group group = file.openGroup(scalarsGroupName);
  writeScalar(group, name, value);
}

/// Read a scalar value from an HDF5 group
/// \tparam T The type of the scalar value
/// \param group The HDF5 group to read from
/// \param name The name of the attribute
/// \param value The scalar value to read
template <typename T>
void readScalar(H5::Group &group, const std::string &name, T &value) {
  H5::Attribute attribute = group.openAttribute(name);
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
  if constexpr (std::is_same_v<T, int>) {
    datatype = H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, double>) {
    datatype = H5::PredType::NATIVE_DOUBLE;
  } // Add more types as needed

  attribute.read(datatype, &value);
}

/// Read a scalar value to an HDF5 file
/// \tparam T The type of the scalar value
/// \param filename The HDF5 file to write to, defoulting to group
/// scalarGroupName
/// \param name The name of the attribute
/// \param value The scalar value to write
template <typename T>
void readScalar(const std::string &filename, const std::string &name,
                T &value) {
  H5::H5File file(filename, H5F_ACC_RDONLY);
  H5::Group group = file.openGroup(scalarsGroupName);
  readScalar(group, name, value);
}

/// @brief write an STL vector to Van HDF5 group
/// @tparam T The type of the vector elements
/// @param group THe HDF5 group to write to
/// @param name  The name of the dataset
/// @param vec  The vector to write
template <typename T>
void writeVector(H5::Group &group, const std::string &name,
                 const std::vector<T> &vec) {
  hsize_t dims[1] = {vec.size()};
  H5::DataSpace dataspace(1, dims);
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
  if constexpr (std::is_same_v<T, int>) {
    datatype = H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, double>) {
    datatype = H5::PredType::NATIVE_DOUBLE;
  } // Add more types as needed

  H5::DataSet dataset = group.createDataSet(name, datatype, dataspace);
  dataset.write(vec.data(), datatype);
}
/// @brief write a C-array to an HDF5 group
/// @tparam T The type of the array elements
/// @param group The HDF5 group to write to
/// @param name The name of the dataset
/// @param vec The array to write
/// @param sz The size of the array
template <typename T>
void writeVector(H5::Group &group, const std::string &name, const T *vec,
                 const int sz) {
  hsize_t dims[1];
  dims[0] = static_cast<hsize_t>(sz);
  H5::DataSpace dataspace(1, dims);
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
  if constexpr (std::is_same_v<T, int>) {
    datatype = H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, double>) {
    datatype = H5::PredType::NATIVE_DOUBLE;
  } // Add more types as needed

  H5::DataSet dataset = group.createDataSet(name, datatype, dataspace);
  dataset.write(vec, datatype);
}

/// @brief read a C-array from an HDF5 group
/// @tparam T The type of the array elements
/// @param group The HDF5 group to write to
/// @param name The name of the dataset
/// @param vec The array to read
/// @param sz The size of the array
template <typename T>
void readVector(H5::Group &group, const std::string &name, T *vec,
                const int sz) {
  H5::DataSet dataset = group.openDataSet(name);
  H5::DataSpace dataspace = dataset.getSpace();
  hsize_t dims[1];
  int ndims = dataspace.getSimpleExtentDims(dims, NULL);

  // if (ndims != 1) {
  //     throw std::runtime_error("Data is not 1-dimensional");
  // }
  assert(ndims == 1);
  assert(dims[0] == sz);

  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
  if constexpr (std::is_same_v<T, int>) {
    datatype = H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, double>) {
    datatype = H5::PredType::NATIVE_DOUBLE;
  } // Add more types as needed

  dataset.read(vec, datatype);
}

/// @brief read an STL vector from an HDF5 group
/// @tparam T The type of the vector elements
/// @param group The HDF5 group to read from
/// @param name The name of the dataset
/// @param vec The vector to read
template <typename T>
void readVector(H5::Group &group, const std::string &name,
                std::vector<T> &vec) {
  H5::DataSet dataset = group.openDataSet(name);
  H5::DataSpace dataspace = dataset.getSpace();
  hsize_t dims[1];
  int ndims = dataspace.getSimpleExtentDims(dims, NULL);

  if (ndims != 1) {
    throw std::runtime_error("Data is not 1-dimensional");
  }

  vec.resize(dims[0]);
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
  if constexpr (std::is_same_v<T, int>) {
    datatype = H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, double>) {
    datatype = H5::PredType::NATIVE_DOUBLE;
  } // Add more types as needed

  dataset.read(vec.data(), datatype);
}


// Helper template to map C++ types to HDF5 datatypes.
template<typename T>
DataType getHDF5Datatype();

template<>
inline DataType getHDF5Datatype<int>() {
    return H5::PredType::NATIVE_INT;
}

template<>
inline DataType getHDF5Datatype<float>() {
    return H5::PredType::NATIVE_FLOAT;
}

template<>
inline DataType getHDF5Datatype<double>() {
    return H5::PredType::NATIVE_DOUBLE;
}

// Add more specializations as needed for other types.
template <typename ViewType>
void readView(ViewType& view, const std::string& fileName, const std::string& name) {
    // Ensure the view resides in host memory.
    static_assert(std::is_same<typename ViewType::memory_space, Kokkos::HostSpace>::value,
                  "The provided view must reside in HostSpace");

    // Open the file for reading.
    H5::H5File file(fileName, H5F_ACC_RDONLY);

    // Open the "views" group.
    H5::Group group;
    try {
        group = file.openGroup(viewsGroupName);
    } catch (H5::Exception &e) {
        throw std::runtime_error("The 'views' group was not found in the file.");
    }

    // Determine the dataset name.
    std::string datasetName = name;
    if (name == "") {
        datasetName = view.label();
    }

    // Open the dataset.
    H5::DataSet dataset = group.openDataSet(datasetName);

    // Retrieve the dataspace and check rank.
    H5::DataSpace dataspace = dataset.getSpace();
    constexpr int rank = ViewType::rank;
    int ndims = dataspace.getSimpleExtentNdims();
    if (ndims != rank) {
        std::cerr << "Name: " << datasetName << "\n";
        throw std::runtime_error("Rank mismatch between dataset and view.");
    }

    // Get dimensions from the dataset.
    hsize_t dims[rank];
    dataspace.getSimpleExtentDims(dims, nullptr);

    // Verify that the dataset dimensions match those of the view.
    bool dimsMatch = true;
    for (int i = 0; i < rank; ++i) {
        if (view.extent(i) != dims[i])
            dimsMatch = false;
    }
    if (!dimsMatch) {
        throw std::runtime_error("Dataset dimensions do not match view dimensions.");
    }

    // Compute the total number of elements.
    size_t total_size = 1;
    for (int i = 0; i < rank; ++i) {
        total_size *= dims[i];
    }

    // Create an std::vector to hold the data.
    typedef typename ViewType::value_type ValueType;
    std::vector<ValueType> buffer(total_size);

    // Create an unmanaged view on the std::vector.
    H5::DataType h5Type = getHDF5Datatype<ValueType>();

    // Read data from the dataset into the unmanaged view.
    dataset.read(buffer.data(), h5Type);

    // Copy data from the unmanaged view into the provided view.
    // Since both views are contiguous and of the same size, a simple copy works.
    for (int i = 0; i < buffer.size(); ++i) {
      int indeces[rank];
      indeces[rank-1] = i;
      for (int ii = rank-1; ii > 0; --ii) {
        indeces[ii-1] = indeces[ii] / dims[ii];
        indeces[ii] %= dims[ii];
      }
      if constexpr (rank == 1) view(indeces[0]) = buffer[i];
      if constexpr (rank == 2) view(indeces[0], indeces[1]) = buffer[i];
      if constexpr (rank == 3) view(indeces[0], indeces[1], indeces[2]) = buffer[i];
      if constexpr (rank == 4) view(indeces[0], indeces[1], indeces[2], indeces[3]) = buffer[i];
      if constexpr (rank == 5) view(indeces[0], indeces[1], indeces[2], indeces[3], indeces[4]) = buffer[i];
      if constexpr (rank == 6) view(indeces[0], indeces[1], indeces[2], indeces[3], indeces[4], indeces[5]) = buffer[i];
      if constexpr (rank == 7) view(indeces[0], indeces[1], indeces[2], indeces[3], indeces[4], indeces[5], indeces[6]) = buffer[i];
    }
}

template <typename ViewType>
void writeView(const ViewType& view, const std::string& fileName, const std::string& name = "") {
    // Ensure the view is in host memory.
    static_assert(std::is_same<typename ViewType::memory_space, Kokkos::HostSpace>::value,
                  "The provided view must reside in HostSpace");

    // Open the file if it exists, or create a new file otherwise.
    H5::H5File file;
    try {
        file = H5::H5File(fileName, H5F_ACC_RDWR);
    } catch (Exception &e) {
        file = H5::H5File(fileName, H5F_ACC_TRUNC);
    }

    // Open the group "views" if it exists; otherwise, create it.
    H5::Group group;
    try {
        group = file.openGroup(viewsGroupName);
    } catch (Exception& e) {
        group = file.createGroup(viewsGroupName);
    }
		std::string datasetName = name;
    // Use the view's label as the dataset name.
    if (name == "") datasetName = view.label();

    // If a dataset with the same name already exists, remove it.
    try {
        H5::DataSet existingDataset = group.openDataSet(datasetName);
        group.unlink(datasetName);
    } catch (H5::Exception &e) {
        // Dataset does not exist; no action needed.
    }


    // Prepare the dataspace using the view dimensions.
    // If the view uses LayoutLeft (column-major), reverse the dimension order.
    constexpr int rank = ViewType::rank;
    hsize_t dims[rank];
    if (std::is_same<typename ViewType::array_layout, Kokkos::LayoutLeft>::value) {
        for (int i = 0; i < rank; ++i) {
            dims[i] = view.extent(rank - i - 1);
        }
    } else {
        // For LayoutRight (row-major) or any other contiguous layout, use the natural order.
        for (int i = 0; i < rank; ++i) {
            dims[i] = view.extent(i);
        }
    }
    H5::DataSpace dataspace(rank, dims);


    // Determine the HDF5 datatype corresponding to the view's value type.
    H5::DataType h5Type = getHDF5Datatype<typename ViewType::value_type>();

    // Create the dataset under the "views" group.
    H5::DataSet dataset = group.createDataSet(datasetName, h5Type, dataspace);

    // Write the data from the view into the dataset.
    dataset.write(view.data(), h5Type);

    // Resources (dataset, group, file) will be closed automatically upon destruction.
}

} // namespace H5
