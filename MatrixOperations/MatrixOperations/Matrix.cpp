#include <iostream>
#include <time.h>
#include <iomanip>
using namespace std;

class Matrix
{
private:

	double** SquareMatrix;

	int LengthOfMatrix, HeightOfMatrix;

	void InitializeSquareMatrix()
	{
		SquareMatrix = new double*[HeightOfMatrix];

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			SquareMatrix[CoordinateX] = new double[LengthOfMatrix];
		}

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix; CoordinateY++)
			{
				SquareMatrix[CoordinateX][CoordinateY] = 0;
			}
		}
	}

public:

	Matrix()
	{
		SquareMatrix = NULL;

		LengthOfMatrix = 0;

		HeightOfMatrix = 0;
	}

	Matrix(int SizeOfMatrix)
	{
		this->LengthOfMatrix = SizeOfMatrix + 1;

		this->HeightOfMatrix = SizeOfMatrix;

		InitializeSquareMatrix();
	}

	int GetLengthOfMatrix()
	{
		return LengthOfMatrix;
	}

	int GetHeightOfMatrix()
	{
		return HeightOfMatrix;
	}

	void DisposeMatrix(const Matrix& Matrix)
	{
		for (int CoordinateX = 0; CoordinateX < Matrix.HeightOfMatrix; CoordinateX++)
		{
			delete[] Matrix.SquareMatrix[CoordinateX];
		}

		delete[] Matrix.SquareMatrix;
	}

	void AsimmetricalMatrixFilling()
	{
		srand(time(0));

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix; CoordinateY++)
			{
				SquareMatrix[CoordinateX][CoordinateY] = 1 + rand() % 100;
			}
		}
	}

	void SimmetricalMatrixFilling()
	{
		srand(time(0));

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix - 1; CoordinateY++)
			{
				if (CoordinateX == CoordinateY)
				{
					SquareMatrix[CoordinateX][CoordinateY] = 1 + rand() % 100;
				}
				else
				{
					SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateY][CoordinateX] = 1 + rand() % 100;
				}
			}

			SquareMatrix[CoordinateX][LengthOfMatrix - 1] = 1 + rand() % 100;
		}
	}

	double* operator [] (int Value)
	{
		return SquareMatrix[Value];
	}

	Matrix operator ~ ()
	{
		double** AxiliaryMatrix = SquareMatrix;

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix; CoordinateY++)
			{
				SquareMatrix[CoordinateY][CoordinateX - 1 - CoordinateX] = AxiliaryMatrix[CoordinateX][CoordinateY];
			}
		}

		return *this;
	}

	Matrix operator * (const double Coefficient)
	{
		Matrix MultiplicationOfMatrix(HeightOfMatrix);

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix; CoordinateY++)
			{
				MultiplicationOfMatrix.SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateX][CoordinateY] * Coefficient;
			}
		}

		return MultiplicationOfMatrix;
	}

	Matrix operator = (const Matrix& SecondMatrix)
	{
		if (SquareMatrix != NULL)
		{
			for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
			{
				delete[] SquareMatrix[CoordinateX];
			}

			delete[] SquareMatrix;
		}

		SquareMatrix = SecondMatrix.SquareMatrix;

		HeightOfMatrix = SecondMatrix.HeightOfMatrix;

		LengthOfMatrix = SecondMatrix.LengthOfMatrix;

		return SecondMatrix;
	}

	Matrix operator * (const Matrix& SecondMatrix)
	{
		if (HeightOfMatrix == SecondMatrix.HeightOfMatrix && LengthOfMatrix == SecondMatrix.LengthOfMatrix)
		{
			Matrix MultilpicationOfMatrix(SecondMatrix.HeightOfMatrix);

			for (int AdditionalCoordinate = 0; AdditionalCoordinate < SecondMatrix.HeightOfMatrix; AdditionalCoordinate++)
			{
				for (int CoordinateX = 0; CoordinateX < SecondMatrix.HeightOfMatrix; CoordinateX++)
				{
					for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
					{
						MultilpicationOfMatrix.SquareMatrix[AdditionalCoordinate][CoordinateX] = SquareMatrix[AdditionalCoordinate][CoordinateY] * SecondMatrix.SquareMatrix[CoordinateY][CoordinateX];
					}
				}
			}

			return MultilpicationOfMatrix;
		}
	}

	Matrix operator - (const Matrix& SecondMatrix)
	{
		if (HeightOfMatrix == SecondMatrix.HeightOfMatrix && LengthOfMatrix == SecondMatrix.LengthOfMatrix)
		{
			Matrix SubtractionOfMatrix(SecondMatrix.HeightOfMatrix);

			for (int CoordinateX = 0; CoordinateX < SecondMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
				{
					SubtractionOfMatrix.SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateX][CoordinateY] - SecondMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}

			return SubtractionOfMatrix;
		}
	}

	Matrix operator + (const Matrix& SecondMatrix)
	{
		if (HeightOfMatrix == SecondMatrix.HeightOfMatrix && LengthOfMatrix == SecondMatrix.LengthOfMatrix)
		{
			Matrix SummitionOfMatrix(SecondMatrix.HeightOfMatrix);

			for (int CoordinateX = 0; CoordinateX < SecondMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
				{
					SummitionOfMatrix.SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateX][CoordinateY] + SecondMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}

			return SummitionOfMatrix;
		}
	}

	friend bool operator == (const Matrix& FirstMatrix, const Matrix& SecondMatrix)
	{
		if (FirstMatrix.HeightOfMatrix != SecondMatrix.HeightOfMatrix || FirstMatrix.LengthOfMatrix != SecondMatrix.LengthOfMatrix)
		{
			return false;
		}

		for (int CoordinateX = 0; CoordinateX < FirstMatrix.HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
			{
				if (FirstMatrix.SquareMatrix[CoordinateX][CoordinateY] != SecondMatrix.SquareMatrix[CoordinateX][CoordinateY])
				{
					return false;
				}
			}
		}
	}

	friend ostream& operator << (ostream& OutputStream, const Matrix& Matrix)
	{
		for (int CoordinateX = 0; CoordinateX < Matrix.HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < Matrix.LengthOfMatrix; CoordinateY++)
			{
				OutputStream << setw(3) << setprecision(3) << Matrix.SquareMatrix[CoordinateX][CoordinateY] << "*x" << CoordinateY;

				if (CoordinateY < Matrix.LengthOfMatrix - 1)
				{
					OutputStream << setw(3) << setprecision(3) << " + ";
				}
				else if(CoordinateY == Matrix.LengthOfMatrix - 1)
				{
					OutputStream << setw(3) << setprecision(3) << "=" <<Matrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}

			OutputStream << endl;
		}

		return OutputStream;
	}

	friend istream& operator >> (istream& InputStream, const double Value)
	{
		InputStream;
	}
};

void GaussMethod(const Matrix& UserMatrix)
{
	Matrix AxillaryMatrix = UserMatrix,ValidationMatrix = UserMatrix;

	int CoordinateX, CoordinateY, AdditionalCoordinate, 

	SizeOfSolutionsArray = AxillaryMatrix.GetHeightOfMatrix();

	double* ArrayOfSolutions = new double[SizeOfSolutionsArray];

	double* ValidationCheck = new double[SizeOfSolutionsArray];

	double Buffer = 0;

	for (CoordinateX = 0; CoordinateX < AxillaryMatrix.GetHeightOfMatrix() - 1; CoordinateX++)
	{
		for (CoordinateY = CoordinateX + 1; CoordinateY < AxillaryMatrix.GetHeightOfMatrix(); CoordinateY++)
		{
			Buffer = AxillaryMatrix[CoordinateX][CoordinateX] / AxillaryMatrix[CoordinateY][CoordinateX];

			for (AdditionalCoordinate = 0; AdditionalCoordinate < AxillaryMatrix.GetHeightOfMatrix(); AdditionalCoordinate++)
			{
				AxillaryMatrix[CoordinateY][AdditionalCoordinate] = AxillaryMatrix[CoordinateY][AdditionalCoordinate] * Buffer - AxillaryMatrix[CoordinateX][AdditionalCoordinate];
			}
		}
	}

	ArrayOfSolutions[SizeOfSolutionsArray - 1] = AxillaryMatrix[SizeOfSolutionsArray - 1][SizeOfSolutionsArray] / AxillaryMatrix[SizeOfSolutionsArray - 1][SizeOfSolutionsArray - 1];

	for (CoordinateX = AxillaryMatrix.GetHeightOfMatrix() - 2; CoordinateX >= 0; CoordinateX--)
	{
		Buffer = 0;

		for (CoordinateY = CoordinateX + 1; CoordinateY < AxillaryMatrix.GetHeightOfMatrix(); CoordinateY++)
		{
			Buffer += AxillaryMatrix[CoordinateX][CoordinateY] * ArrayOfSolutions[CoordinateY];
		}

		ArrayOfSolutions[CoordinateX] = (AxillaryMatrix[CoordinateX][SizeOfSolutionsArray] - Buffer) / AxillaryMatrix[CoordinateX][CoordinateX];
	}

	cout << "A solution of current system of linear equations :" << endl;

	for (CoordinateX = 0; CoordinateX < SizeOfSolutionsArray; CoordinateX++)
	{
		cout << "x" << CoordinateX << " = " << ArrayOfSolutions[CoordinateX] << endl;
	}

	cout << "Validation checking :" << endl;

	for (CoordinateX = 0; CoordinateX < SizeOfSolutionsArray; CoordinateX++)
	{
		ValidationCheck[CoordinateX] = 0;

		for (CoordinateY = 0; CoordinateY < SizeOfSolutionsArray; CoordinateY++)
		{
			ValidationCheck[CoordinateX] += ValidationMatrix[CoordinateX][CoordinateY] * ArrayOfSolutions[CoordinateY];
		}

		cout << ValidationCheck[CoordinateX] << endl;
	}

	delete[] ValidationCheck;

	delete[] ArrayOfSolutions;
}

void main()
{
	Matrix UserMatrix(7);

	UserMatrix.SimmetricalMatrixFilling();

	cout << UserMatrix;

	GaussMethod(UserMatrix);

	system("pause");
}