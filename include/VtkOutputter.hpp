// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file VtkOutputter.hpp
 * @brief Defines the VtkOutputter class that writes output to *.vtu files.
 */
#pragma once

#ifndef __MINGW32__

#include <vtk-7.1/vtkXMLUnstructuredGridWriter.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <memory>
#include "Types.hpp"
#include "TriangleNodesProvisioner.hpp"
#include "OutputterBase.hpp"
#include <boost/python/numpy.hpp>
#include <boost/python.hpp>

namespace blitzdg {
  /**
   * Outputter class for vtk files.
   */ 
  class VtkOutputter : OutputterBase {

	vtkSmartPointer<vtkXMLUnstructuredGridWriter> GridWriter;
	const TriangleNodesProvisioner& NodesProvisioner;
	std::string FileExtension;

public:
	VtkOutputter(TriangleNodesProvisioner & _NodesProvisioner);

	std::string generateFileName(const std::string & fieldName, const index_type fileNumber) const;

	void writeFieldToFile(const std::string & fileName, real_matrix_type field, const std::string & fieldName) const;
	void writeFieldsToFiles(std::map<std::string, real_matrix_type>& fields, index_type tstep);

	void writeFieldToFile_numpy(boost::python::str fileName, const boost::python::numpy::ndarray& field, boost::python::str fieldName) const;
	void writeFieldsToFiles_numpy(const boost::python::dict& fields, index_type tstep);



private:
	void splitTriangles(const real_matrix_type& x, const real_matrix_type& y, const real_matrix_type& field, real_matrix_type& xnew, real_matrix_type& ynew, real_matrix_type& fieldnew) const;

  };
}
#endif
