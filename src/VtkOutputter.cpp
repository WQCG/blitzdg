// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#ifndef __MINGW32__

#include "Types.hpp"
#include "VtkOutputter.hpp"
#include "BlitzHelpers.hpp"
#include <vtk-7.1/vtkXMLUnstructuredGridWriter.h>
#include <vtk-7.1/vtkIndent.h>
#include <vtk-7.1/vtkUnstructuredGrid.h>
#include <vtk-7.1/vtkPoints.h>
#include <vtk-7.1/vtkProperty.h>
#include <vtk-7.1/vtkType.h>
#include <vtk-7.1/vtkTriangle.h>
#include <vtk-7.1/vtkDataSetAttributes.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkCellArray.h>
#include <vtk-7.1/vtkDoubleArray.h>
#include <vtk-7.1/vtkFieldData.h>
#include <vtk-7.1/vtkPointData.h>
#include <vtk-7.1/vtkPolyData.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <boost/python.hpp>

using std::stringstream;
using std::ostream;
using std::ios;
using std::filebuf;
using std::string;
using std::vector;
using boost::python::str;
using boost::python::numpy::ndarray;
using boost::python::stl_input_iterator;

namespace blitzdg {
	VtkOutputter::VtkOutputter(NodesProvisioner2DBase& _NodesProvisioner)
		: GridWriter { vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New() },
			NodesProvisioner { _NodesProvisioner } 
	{
		FileExtension = GridWriter->GetDefaultFileExtension();
	}

	/**
	 * Generates a file name for storing a blitzdg object in vtk format.
	 * @param[in] fieldName Name of field from the PDE (system), e.g., "u".
	 * @param[in] fileNumber An integral intex indicating a logical ordering on the output files. It is usually related to time-level.
	 */
	string VtkOutputter::generateFileName(const string & fieldName, const index_type fileNumber) const {
		stringstream fileNameStrm;
		fileNameStrm << fieldName << setfill('0') << setw(7) << fileNumber << "." << FileExtension;
		return fileNameStrm.str();
	}

	void VtkOutputter::writeFieldsToFiles(std::map<std::string, real_matrix_type>& fields, index_type tstep) {
		for (const auto& kv : fields) {
			const string& fieldName = kv.first;
			const real_matrix_type& field = kv.second;
			string fileName = generateFileName(fieldName, tstep);
			writeFieldToFile(fileName, field, fieldName);
		}
	}


	void VtkOutputter::writeFieldToFile_numpy(boost::python::str fileName, const boost::python::numpy::ndarray& field, boost::python::str fieldName) const
 	{
		real_matrix_type blitzField(field.shape(0), field.shape(1));
		blitzField = 0.0;
		char * raw = field.get_data();
        std::copy(&raw[0], &raw[field.shape(0)*field.shape(1)*sizeof(real_type)] , reinterpret_cast<char*>(blitzField.data()));

		string fileNameCpp = string(boost::python::extract<const char*>(fileName));
		string fieldNameCpp = string(boost::python::extract<const char*>(fieldName));

		writeFieldToFile(fileNameCpp, blitzField, fieldNameCpp);
	}

	void VtkOutputter::writeFieldsToFiles_numpy(const boost::python::dict& fields, index_type tstep) {

		auto fieldNames = std::vector<const char*>(boost::python::stl_input_iterator<const char *>(fields.keys()),
			stl_input_iterator<const char *>() );

		auto fieldArrays = std::vector<ndarray>(boost::python::stl_input_iterator<ndarray>(fields.values()),
			stl_input_iterator<ndarray>()); 

		index_type count = 0;
		for (const auto& f : fieldNames) {
			const char * fieldNameC = f;
			string fieldNameCpp = fieldNameC;
			string fileName = generateFileName(fieldNameCpp, tstep);
			const ndarray& field = fieldArrays.at(count);
			
			writeFieldToFile_numpy(str(fileName), field, str(fieldNameCpp));
			++count;
		}
	}
}
#endif