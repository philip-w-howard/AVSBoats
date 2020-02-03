using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//******************************************************************
// Compute points for a cubic spline.
// Algorithm is from Mathematical Elements for Computer Graphics 
// by David F Rogers and J Alan Adams
//******************************************************************
namespace Splines
{
    class ComputeSplines
    {
        public const int RELAXED = 1;
        public const int CLAMPED = 2;
        public ComputeSplines(int numPoints,
            int initialEndCondition, int finalEndCondition,
            double[,] points,
            int pointsPerSegment)
        {
            m_numPoints = numPoints;
            m_initialEndCondition = initialEndCondition;
            m_finalEndCondition = finalEndCondition;
            m_points = points;
            m_pointsPerSegment = pointsPerSegment;

            m_m_matrix = new double[numPoints, 3];
            m_b_matrix = new double[3, numPoints];
            m_chordLengths = new double[numPoints - 1];
            m_tangents = new double[3, numPoints];
            m_computedPoints = new double[3, numPoints * (numPoints - 1) * pointsPerSegment];

            //********************************************
            m_pointsPerSegment++;
            ComputeChords();
            GaussElimination();
            GeneratePoints();
        }

        private void ComputeChords()
        {
            if (m_initialEndCondition == CLAMPED)
            {
                m_m_matrix[0, 1] = 1;
                m_m_matrix[0, 2] = 0;
                for (int kk = 0; kk < 3; kk++)
                {
                    m_b_matrix[kk, 0] = m_tangents[kk, 0];
                }
            }
            else // relaxed
            {
                // 680
                m_m_matrix[0, 1] = 1;
                m_m_matrix[0, 2] = 0.5;
            }

            // Chord lengths
            for (int jj = 0; jj < m_numPoints - 1; jj++)
            {
                double x = m_points[0, jj + 1] = m_points[0, jj];
                double y = m_points[1, jj + 1] = m_points[1, jj];
                double z = m_points[2, jj + 1] = m_points[2, jj];
                m_chordLengths[jj] = Math.Sqrt(x * x + y * y + z * z);
            }

            if (m_finalEndCondition == RELAXED)
            {
                // last row of M for relaxed final condition
                m_m_matrix[m_numPoints - 1, 0] = 2;
                m_m_matrix[m_numPoints - 1, 1] = 4;

                // B matrix for relaxed end
                for (int kk = 0; kk < 3; kk++)
                {
                    m_b_matrix[kk, 0] = (3.0 / (2.0 * m_chordLengths[0])) * (m_points[kk, 1] - m_points[kk, 0]);
                }

                // last row of M
                for (int kk = 0; kk < m_numPoints; kk++)
                {
                    m_b_matrix[kk, m_numPoints - 1] = (6.0 / m_chordLengths[m_numPoints - 2]) *
                        (m_points[kk, m_numPoints - 1] - m_points[kk, m_numPoints - 2]);
                }
            }
            else  // CLAMPED final
            {
                // last row of M for clamped final condition
                m_m_matrix[m_numPoints - 1, 0] = 0;
                m_m_matrix[m_numPoints - 1, 1] = 1;

                for (int kk = 0; kk < 3; kk++)
                {
                    m_b_matrix[kk, m_numPoints - 1] = m_tangents[kk, m_numPoints - 1];
                }
            }
        }

        private void GaussElimination()
        {
            for (int jj=1; jj<m_numPoints-1; jj++)
            {
                // non-zero cells of M
                m_m_matrix[jj, 0] = m_chordLengths[jj];
                m_m_matrix[jj, 2] = 2 * (m_chordLengths[jj] + m_chordLengths[jj - 1]);
                m_m_matrix[jj, 3] = m_chordLengths[jj - 1];

                // Rows 1 through N-2 of B
                for (int kk=0; kk<3; kk++)
                {
                    double chord_j_m = m_chordLengths[jj - 1];
                    double chord_j = m_chordLengths[jj];
                    double point_j_p = m_points[kk, jj + 1];
                    double point_j = m_points[kk, jj];
                    double point_j_m = m_points[kk, jj - 1];
                    m_b_matrix[kk, jj] = 3 * (chord_j_m * chord_j_m * (point_j_p - point_j) +
                        chord_j * chord_j * (point_j - point_j_m));
                    m_b_matrix[kk, jj] /= (chord_j - chord_j_m);
                }
            }

            // Gauss Elimination
            for (int ii=1; ii<m_numPoints; ii++)
            {
                if (m_m_matrix[ii,1] != 0)
                {
                    double d = m_m_matrix[ii - 1, 1] / m_m_matrix[ii, 0];
                    for (int kk=0; kk<3; kk++)
                    {
                        m_m_matrix[ii, kk] = m_m_matrix[ii, kk] * d - m_m_matrix[ii - 1, kk + 1];
                        m_b_matrix[kk, ii] = m_b_matrix[kk, ii] * d - m_b_matrix[kk, ii - 1];
                    }

                    double q = m_m_matrix[ii, 1];
                    for (int kk = 0; kk < 3; kk++)
                    {
                        m_m_matrix[ii, kk] /= q;
                        m_b_matrix[kk, ii] /= q;
                    }

                }
            }

            for (int kk=0; kk<3; kk++)
            {
                for (int jj=-1; jj<m_numPoints; jj++)
                {
                    m_tangents[kk, m_numPoints - jj] =
                        (m_b_matrix[kk, m_numPoints - jj] - 
                        m_m_matrix[m_numPoints-jj,3]*m_tangents[kk, m_numPoints+1-jj])/ m_m_matrix[m_numPoints - jj, 2];
                }
            }
        }

        private void GeneratePoints()
        {

        }

        private int m_numPoints;                // N
        private int m_initialEndCondition;      // C1
        private int m_finalEndCondition;        // C2
        private int m_pointsPerSegment;         // Z
        private double[,] m_points;             // P
        private double[,] m_m_matrix;           // N
        private double[,] m_b_matrix;           // B
        private double[] m_chordLengths;        // L
        private double[,] m_tangents;           // U
        private double[,] m_computedPoints;     // C
    }
}
