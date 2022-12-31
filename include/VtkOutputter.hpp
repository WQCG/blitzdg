// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file VtkOutputter.hpp
 * @brief Defines the VtkOutputter class that writes output to *.vtu files.
 */
#pragma once

#ifndef __MINGW32__

#include <vtk-7.1/vtkXMLUnstructuredGridWriter.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkPointData.h>
#include <vtk-7.1/vtkCellArray.h>
#include <vtk-7.1/vtkDoubleArray.h>
#include <vtk-7.1/vtkUnstructuredGrid.h>
#include <boost/python/numpy.hpp>
#include <boost/python.hpp>

#include <memory>
#include "Types.hpp"
#include "TriangleNodesProvisioner.hpp"
#include "OutputterBase.hpp"


namespace blitzdg {
  /**
   * Outputter class for vtk files.
   */ 
  class VtkOutputter : OutputterBase {

	vtkSmartPointer<vtkXMLUnstructuredGridWriter> GridWriter;
	const NodesProvisioner2DBase& NodesProvisioner;
	std::string FileExtension;

public:
	VtkOutputter(NodesProvisioner2DBase & _NodesProvisioner);

	std::string generateFileName(const std::string & fieldName, const index_type fileNumber) const;

	/**
     * Writes a blitz array to plain-text file.
     * @param[in] fileName Name of the file (e.g., field0000010.vtu).
     * @param[in] field Two-dimensional blitz array to be written to the file. Usually a 'field' of the PDE (system) being solved.
     * @param[in] delimeter Character that will be used to separate columns. Rows are always separated by line-endings.
     */
    void writeFieldToFile(const std::string & fileName, real_matrix_type field, const std::string & fieldName) const {

		real_matrix_type x = NodesProvisioner.get_xGrid(), y = NodesProvisioner.get_yGrid();

		index_type K = field.cols();
		index_type Np = field.rows();
		// If higher order than linear, need to break up the trangles.
		if (NodesProvisioner.get_NOrder() > 1) {
			real_matrix_type xnew, ynew, fieldnew;
			NodesProvisioner.splitElements(x, y, field, xnew, ynew, fieldnew);

			Np = fieldnew.rows();
			K = fieldnew.cols();

			x.resize(Np, K);
			y.resize(Np, K);
			field.resize(Np, K);

			x = xnew;
			y = ynew;
			field = fieldnew;
		}

		vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

		vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();

		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

		const char* fieldNameChar = fieldName.c_str();
		array->SetName(fieldNameChar);
		array->SetNumberOfValues(Np*K);

		unstructuredGrid->Allocate(K);		
		insertAllCells(x, y, field, points, array, unstructuredGrid);

		unstructuredGrid->SetPoints(points);
		unstructuredGrid->GetPointData()->SetScalars(array);
		unstructuredGrid->GetPointData()->SetActiveScalars(fieldNameChar);

		// Write file
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName(fileName.c_str());
#if VTK_MAJOR_VERSION <= 5
		writer->SetInput(unstructuredGrid);
#else
		writer->SetInputData(unstructuredGrid);
#endif
		writer->Write();
    }
	
	void writeFieldsToFiles(std::map<std::string, real_matrix_type>& fields, index_type tstep);

	void writeFieldToFile_numpy(boost::python::str fileName, const boost::python::numpy::ndarray& field, boost::python::str fieldName) const;
	void writeFieldsToFiles_numpy(const boost::python::dict& fields, index_type tstep);



private:
	void splitTriangles(const real_matrix_type& x, const real_matrix_type& y, const real_matrix_type& field, real_matrix_type& xnew, real_matrix_type& ynew, real_matrix_type& fieldnew) const;

	void insertAllCells(real_matrix_type& x, real_matrix_type& y, real_matrix_type& field, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> array, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) const {
            // '3' because Triangles
            // vtkIdType nodes[3];
            index_type K  = field.cols();
            index_type Np = field.rows();
            std::vector<vtkIdType> nodes(Np), permutedNodes(Np);
            std::map<index_type, VTKCellType> cellTypes = 
                {
                    {3, VTK_TRIANGLE},
                    {4, VTK_QUAD}
                };

            index_type nodeId = 0;
            for (index_type k=0; k < K; ++k) {
                for (index_type n=0; n < Np; ++n) {
                    points->InsertPoint(nodeId, x(n, k), y(n, k), 0.0 );
                    
                    array->SetValue(nodeId, field(n,k));

                    nodes[n] = nodeId;
                    ++nodeId;
                }

                if (cellTypes[Np] == VTK_QUAD) {
					std::swap(nodes[1], nodes[3]);
					std::swap(nodes[1], nodes[2]);
                }

				unstructuredGrid->InsertNextCell(cellTypes[Np], Np, nodes.data());
            }
        }

  };
}
#endif
