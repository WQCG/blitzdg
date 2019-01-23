#ifndef __MINGW32__

#include "Types.hpp"
#include "VtkOutputter.hpp"
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

using std::stringstream;
using std::ostream;
using std::ios;
using std::filebuf;
using std::string;

namespace blitzdg{
	VtkOutputter::VtkOutputter(TriangleNodesProvisioner & _NodesProvisioner)
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
		string VtkOutputter::generateFileName(const string & fieldName, const index_type fileNumber) {
			stringstream fileNameStrm;
        	fileNameStrm << fieldName << setfill('0') << setw(7) << fileNumber << "." << FileExtension;
        	return fileNameStrm.str();
		}


	/**
     * Writes a blitz array to plain-text file.
     * @param[in] fileName Name of the file (e.g., field0000010.vtu).
     * @param[in] field Two-dimensional blitz array to be written to the file. Usually a 'field' of the PDE (system) being solved.
     * @param[in] delimeter Character that will be used to separate columns. Rows are always separated by line-endings.
     */
    void VtkOutputter::writeFieldToFile(const string & fileName, const real_matrix_type & field, const string & fieldName) {

		const real_matrix_type& x = NodesProvisioner.get_xGrid(), y = NodesProvisioner.get_yGrid();

		index_type K = field.cols();
		index_type Np = field.rows();
		index_type nodeId = 0;

		vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

		vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();

		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

		const char* fieldNameChar = fieldName.c_str();
		array->SetName(fieldNameChar);
		array->SetNumberOfValues(Np*K);

		unstructuredGrid->Allocate(K);

		// '3' because this will only work for linear elements until we interpolate
		// to a finer (uniform) triangular mesh at higher order.
		vtkIdType nodes[3];

		for (index_type k=0; k < K; ++k) {
			for (index_type n=0; n < Np; ++n) {
				points->InsertPoint(nodeId, x(n, k), y(n, k), field(n,k) );
				
				array->SetValue(nodeId, field(n,k));

				nodes[n] = nodeId;
				++nodeId;
			}

			unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, nodes);
		}

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
}
#endif