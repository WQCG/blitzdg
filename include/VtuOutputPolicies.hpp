#pragma once

#include <vtk-7.1/vtkUnstructuredGrid.h>
#include <vtk-7.1/vtkXMLUnstructuredGridWriter.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkDoubleArray.h>
#include <vtk-7.1/vtkPoints.h>
#include <vtk-7.1/vtkCellType.h>

#include "Types.hpp"

namespace blitzdg {
    struct TriangleVtuOutputPolicy {
        static void insertAllCells(real_matrix_type& x, real_matrix_type& y, real_matrix_type& field, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> array, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
            // '3' because Triangles
            vtkIdType nodes[3];
            index_type K  = field.cols();
            index_type Np = field.rows();
            index_type nodeId = 0;
            for (index_type k=0; k < K; ++k) {
                for (index_type n=0; n < Np; ++n) {
                    points->InsertPoint(nodeId, x(n, k), y(n, k), field(n,k) );
                    
                    array->SetValue(nodeId, field(n,k));

                    nodes[n] = nodeId;
                    ++nodeId;
                }

                unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, nodes);
            }
        }
    };

    struct QuadVtuOutputPolicy {
        static void insertAllCells(real_matrix_type& x, real_matrix_type& y, real_matrix_type& field, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> array, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid) {
            // '4' because Quadss
            vtkIdType nodes[4];
            index_type K  = field.cols();
            index_type Np = field.rows();
            index_type nodeId = 0;
            for (index_type k=0; k < K; ++k) {
                for (index_type n=0; n < Np; ++n) {
                    points->InsertPoint(nodeId, x(n, k), y(n, k), field(n,k) );
                    
                    array->SetValue(nodeId, field(n,k));

                    nodes[n] = nodeId;
                    ++nodeId;
                }

                unstructuredGrid->InsertNextCell(VTK_QUAD, 4, nodes);
            }
        }
    };
}
