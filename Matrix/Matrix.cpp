#define _USE_MATH_DEFINES
#include <iostream>
#include <time.h>
#include <iomanip>
#include <fstream>
#include <vector>

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

		NullMatrixFilling();
	}

public:

	Matrix()
	{
		SquareMatrix = NULL;

		LengthOfMatrix = 0;

		HeightOfMatrix = 0;
	}

	Matrix(int SizeOfSquareMatrix)
	{
		this->LengthOfMatrix = SizeOfSquareMatrix;

		this->HeightOfMatrix = SizeOfSquareMatrix;

		InitializeSquareMatrix();
	}

	Matrix(int Height, int Length)
	{
		this->HeightOfMatrix = Height;

		this->LengthOfMatrix = Length;

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

	void UnitMatrixFilling()
	{
		int Coefficient = 1;

		if (HeightOfMatrix == LengthOfMatrix)
		{
			Coefficient = 0;
		}

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix - Coefficient; CoordinateY++)
			{
				SquareMatrix[CoordinateX][CoordinateY] = 0;
			}

			SquareMatrix[CoordinateX][CoordinateX] = 1;
		}
	}

	void UnitVectorFilling()
	{
		for (int CoordinateX = 0; CoordinateX < LengthOfMatrix; CoordinateX++)
		{
			SquareMatrix[0][CoordinateX] = 1;
		}
	}

	void NullMatrixFilling()
	{
		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix; CoordinateY++)
			{
				SquareMatrix[CoordinateX][CoordinateY] = 0;
			}
		}
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

		int Coefficient = 1;

		if (HeightOfMatrix == LengthOfMatrix)
		{
			Coefficient = 0;
		}

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix - Coefficient; CoordinateY++)
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

			SquareMatrix[CoordinateX][LengthOfMatrix - Coefficient] = 1 + rand() % 100;
		}
	}

	void HilbertMatrixFilling()
	{
		for (int CoordinateX = 1; CoordinateX <= HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 1; CoordinateY <= LengthOfMatrix; CoordinateY++)
			{
				SquareMatrix[CoordinateX - 1][CoordinateY - 1] = 1.0 / double(CoordinateX + CoordinateY - 1);
			}
		}
	}

	Matrix LinearRegressionMatrixFilling(int NumberOfParameters, Matrix& VectorOfCoefficients)
	{
		Matrix TestMatrix(NumberOfParameters,VectorOfCoefficients.LengthOfMatrix);

		TestMatrix.NullMatrixFilling();

		for (int CoordinateX = 0; CoordinateX < NumberOfParameters; CoordinateX++)
		{
			for (int CoordinateY = 1; CoordinateY < VectorOfCoefficients.LengthOfMatrix; CoordinateY++)
			{
				TestMatrix[CoordinateX][CoordinateY] = (1 + rand() % 1000) * pow(-1,rand());

				TestMatrix[CoordinateX][0] += TestMatrix[CoordinateX][CoordinateY] * VectorOfCoefficients[0][CoordinateY];
			}

			TestMatrix[CoordinateX][0] += VectorOfCoefficients[0][0];

			TestMatrix[CoordinateX][0] += (rand() % 7) * pow(-1, rand());
		}

		return TestMatrix;
	}

	void Dispose()
	{
		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			delete[] SquareMatrix[CoordinateX];
		}

		delete[] SquareMatrix;
	}

	double* operator [] (int Index)
	{
		return SquareMatrix[Index];
	}

	Matrix operator * (double Multiplicator)
	{
		Matrix MultiplicationMatrix(HeightOfMatrix, LengthOfMatrix);

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix; CoordinateY++)
			{
				MultiplicationMatrix.SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateX][CoordinateY] * Multiplicator;
			}
		}

		return MultiplicationMatrix;
	}

	Matrix operator = (const Matrix& SecondMatrix)
	{
		if (HeightOfMatrix == SecondMatrix.HeightOfMatrix && LengthOfMatrix == SecondMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < SecondMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
				{
					SquareMatrix[CoordinateX][CoordinateY] = SecondMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}

		}
		return *this;
	}

	Matrix operator * (Matrix& SecondMatrix)
	{
		Matrix MultiplicationMatrix(SecondMatrix.GetHeightOfMatrix(), SecondMatrix.GetLengthOfMatrix());

		if (HeightOfMatrix == SecondMatrix.HeightOfMatrix && LengthOfMatrix == SecondMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < SecondMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
				{
					for (int AdditionalCoordinate = 0; AdditionalCoordinate < SecondMatrix.LengthOfMatrix; AdditionalCoordinate++)
					{
						MultiplicationMatrix.SquareMatrix[CoordinateX][CoordinateY] += SquareMatrix[CoordinateX][AdditionalCoordinate] * SecondMatrix.SquareMatrix[AdditionalCoordinate][CoordinateY];
					}
				}
			}
		}

		return MultiplicationMatrix;
	}

	Matrix operator - (Matrix& SecondMatrix)
	{
		Matrix SubstractionMatrix(SecondMatrix.GetHeightOfMatrix(), SecondMatrix.GetLengthOfMatrix());

		if (HeightOfMatrix == SecondMatrix.HeightOfMatrix && LengthOfMatrix == SecondMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < SecondMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
				{
					SubstractionMatrix.SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateX][CoordinateY] - SecondMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}

			return SubstractionMatrix;
		}
	}

	Matrix operator + (Matrix& SecondMatrix)
	{
		Matrix SummitionMatrix(SecondMatrix.GetHeightOfMatrix(), SecondMatrix.GetLengthOfMatrix());

		if (HeightOfMatrix == SecondMatrix.HeightOfMatrix && LengthOfMatrix == SecondMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < SecondMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < SecondMatrix.LengthOfMatrix; CoordinateY++)
				{
					SummitionMatrix.SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateX][CoordinateY] + SecondMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}

			return SummitionMatrix;
		}
	}

	friend Matrix operator ~ (Matrix& UserMatrix)
	{
		Matrix TranspositionMatrix(UserMatrix.HeightOfMatrix, UserMatrix.LengthOfMatrix);

		TranspositionMatrix = UserMatrix;

		double Buffer = 0;

		for (int CoordinateX = 0; CoordinateX < UserMatrix.HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < UserMatrix.LengthOfMatrix; CoordinateY++)
			{
				Buffer = TranspositionMatrix.SquareMatrix[CoordinateX][CoordinateY];

				TranspositionMatrix.SquareMatrix[CoordinateX][CoordinateY] = TranspositionMatrix.SquareMatrix[CoordinateY][CoordinateX];

				TranspositionMatrix.SquareMatrix[CoordinateY][CoordinateX] = Buffer;
			}
		}

		return TranspositionMatrix;
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

	friend ostream& operator << (ostream& OutputStream, const Matrix& UserMatrix)
	{
		if (UserMatrix.HeightOfMatrix == UserMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < UserMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < UserMatrix.LengthOfMatrix; CoordinateY++)
				{
					OutputStream << setw(3) << UserMatrix.SquareMatrix[CoordinateX][CoordinateY] << "*x" << CoordinateY + 1;

					if (CoordinateY < UserMatrix.LengthOfMatrix - 1)
					{
						OutputStream << setw(3) << " + ";
					}
				}

				OutputStream << endl;
			}
		}
		else if (UserMatrix.HeightOfMatrix == 1 && UserMatrix.LengthOfMatrix >= 1)
		{
			OutputStream << setw(3) << "x = (";

			for (int CoordinateX = 0; CoordinateX < UserMatrix.LengthOfMatrix; CoordinateX++)
			{
				OutputStream << UserMatrix.SquareMatrix[0][CoordinateX];

				if (CoordinateX < UserMatrix.LengthOfMatrix - 1)
				{
					OutputStream << ",";
				}
			}

			OutputStream << ")" << endl;
		}
		else if (UserMatrix.HeightOfMatrix != UserMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < UserMatrix.HeightOfMatrix; CoordinateX++)
			{
				for (int CoordinateY = 0; CoordinateY < UserMatrix.LengthOfMatrix - 1; CoordinateY++)
				{
					OutputStream << setw(3) << UserMatrix.SquareMatrix[CoordinateX][CoordinateY] << "*x" << CoordinateY + 1;

					if (CoordinateY < UserMatrix.LengthOfMatrix - 2)
					{
						OutputStream << setw(3) << " + ";
					}
				}

				OutputStream << setw(3) << "= " << UserMatrix.SquareMatrix[CoordinateX][UserMatrix.LengthOfMatrix - 1];

				OutputStream << endl;
			}
		}

		return OutputStream;
	}

	friend istream& operator >> (istream& InputStream, const Matrix& UserMatrix)
	{
		if (UserMatrix.HeightOfMatrix == UserMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < UserMatrix.HeightOfMatrix; CoordinateX++)
			{
				cout << "Fill the " << CoordinateX + 1 << " row : ";

				for (int CoordinateY = 0; CoordinateY < UserMatrix.LengthOfMatrix; CoordinateY++)
				{
					InputStream >> UserMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}
		}
		else if (UserMatrix.HeightOfMatrix == 1 && UserMatrix.LengthOfMatrix >= 1)
		{
			for (int CoordinateX = 0; CoordinateX < UserMatrix.HeightOfMatrix; CoordinateX++)
			{
				cout << "Fill the coordinates of vector : ";

				for (int CoordinateY = 0; CoordinateY < UserMatrix.LengthOfMatrix; CoordinateY++)
				{
					InputStream >> UserMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}
		}
		else if (UserMatrix.HeightOfMatrix != UserMatrix.LengthOfMatrix)
		{
			for (int CoordinateX = 0; CoordinateX < UserMatrix.HeightOfMatrix; CoordinateX++)
			{
				cout << "Fill the " << CoordinateX + 1 << " row and " << UserMatrix.LengthOfMatrix << " element has to be free : ";

				for (int CoordinateY = 0; CoordinateY < UserMatrix.LengthOfMatrix; CoordinateY++)
				{
					InputStream >> UserMatrix.SquareMatrix[CoordinateX][CoordinateY];
				}
			}
		}
		cout << endl;

		return InputStream;
	}
};

bool MatrixIsSimmetric(Matrix UserMatrix)
{
	for (int CoordinateX = 0; CoordinateX < UserMatrix.GetHeightOfMatrix(); CoordinateX++)
	{
		for (int CoordinateY = 0; CoordinateY < UserMatrix.GetLengthOfMatrix() - 1; CoordinateY++)
		{
			if (UserMatrix[CoordinateX][CoordinateY] != UserMatrix[CoordinateY][CoordinateX])
			{
				return false;
			}
		}
	}

	return true;
}

double VectorRate(Matrix UserMatrix)
{
	double Rate = 0;

	for (int CoordinateX = 0; CoordinateX < UserMatrix.GetLengthOfMatrix(); CoordinateX++)
	{
		Rate += UserMatrix[0][CoordinateX] * UserMatrix[0][CoordinateX];
	}

	return sqrt(Rate);
}

double ScalarProduct(Matrix FirstVector, Matrix SecondVector)
{
	double Scalar = 0;

	for (int CoordinateX = 0; CoordinateX < FirstVector.GetLengthOfMatrix(); CoordinateX++)
	{
		Scalar += FirstVector[0][CoordinateX] * SecondVector[0][CoordinateX];
	}

	return Scalar;
}

void GaussMethod(Matrix UserMatrix)
{
	Matrix CopiedMatrix = UserMatrix,

		ValidationMatrix(UserMatrix.GetHeightOfMatrix(), UserMatrix.GetLengthOfMatrix());

	int CoordinateX = 0,

		CoordinateY = 0,

		AdditionalCoordinate = 0,

		SizeOfSolutionsArray = CopiedMatrix.GetHeightOfMatrix();

	double* ArrayOfSolutions = new double[SizeOfSolutionsArray];

	double* ValidationArray = new double[SizeOfSolutionsArray];

	double Buffer = 0,

		Result = 0;

	ValidationMatrix = CopiedMatrix;

	cout << "System of linear equations :\n" << endl;

	cout << CopiedMatrix << endl;

	for (CoordinateX = 0; CoordinateX < CopiedMatrix.GetHeightOfMatrix() - 1; CoordinateX++)
	{
		for (CoordinateY = CoordinateX + 1; CoordinateY < CopiedMatrix.GetHeightOfMatrix(); CoordinateY++)
		{
			Buffer = CopiedMatrix[CoordinateX][CoordinateX] / CopiedMatrix[CoordinateY][CoordinateX];

			for (AdditionalCoordinate = 0; AdditionalCoordinate <= CopiedMatrix.GetHeightOfMatrix(); AdditionalCoordinate++)
			{
				CopiedMatrix[CoordinateY][AdditionalCoordinate] = CopiedMatrix[CoordinateY][AdditionalCoordinate] * Buffer - CopiedMatrix[CoordinateX][AdditionalCoordinate];
			}
		}
	}

	ArrayOfSolutions[SizeOfSolutionsArray - 1] = CopiedMatrix[SizeOfSolutionsArray - 1][SizeOfSolutionsArray] / CopiedMatrix[SizeOfSolutionsArray - 1][SizeOfSolutionsArray - 1];

	for (CoordinateX = CopiedMatrix.GetHeightOfMatrix() - 2; CoordinateX >= 0; CoordinateX--)
	{
		Buffer = 0;

		for (CoordinateY = CoordinateX + 1; CoordinateY < CopiedMatrix.GetHeightOfMatrix(); CoordinateY++)
		{
			Buffer += CopiedMatrix[CoordinateX][CoordinateY] * ArrayOfSolutions[CoordinateY];
		}

		ArrayOfSolutions[CoordinateX] = (CopiedMatrix[CoordinateX][SizeOfSolutionsArray] - Buffer) / CopiedMatrix[CoordinateX][CoordinateX];
	}

	cout << "A solution of current system of linear equations :\n" << endl;

	for (CoordinateX = 0; CoordinateX < SizeOfSolutionsArray; CoordinateX++)
	{
		cout << setw(2) << "x" << CoordinateX + 1 << " = " << ArrayOfSolutions[CoordinateX] << endl;

		if (CoordinateX == SizeOfSolutionsArray - 1)
		{
			cout << endl;
		}
	}

	cout << "Validation checking :\n" << endl;

	for (CoordinateX = 0; CoordinateX < SizeOfSolutionsArray; CoordinateX++)
	{
		Result = 0;

		for (CoordinateY = 0; CoordinateY < SizeOfSolutionsArray; CoordinateY++)
		{
			cout << setw(3) << ValidationMatrix[CoordinateX][CoordinateY];

			if (CoordinateY < SizeOfSolutionsArray - 1)
			{
				cout << "*" << ArrayOfSolutions[CoordinateY] << setw(3) << " + ";
			}

			if (CoordinateY == SizeOfSolutionsArray - 1)
			{
				cout << "*" << ArrayOfSolutions[CoordinateY] << setw(3) << " = ";
			}

			Result += ValidationMatrix[CoordinateX][CoordinateY] * ArrayOfSolutions[CoordinateY];
		}

		cout << Result << endl;
	}

	cout << endl;

	delete[] ValidationArray;

	delete[] ArrayOfSolutions;
}

void KaczmarzMethod(Matrix UserMatrix, const double PrecisionOfResult)
{
	Matrix Vector(1, UserMatrix.GetLengthOfMatrix() - 1);

	Matrix AxillaryVector(1, UserMatrix.GetLengthOfMatrix() - 1);

	Matrix SubstractionVector(1, UserMatrix.GetLengthOfMatrix() - 1);

	for (int CoordinateX = 0; CoordinateX < Vector.GetLengthOfMatrix(); CoordinateX++)
	{
		Vector[0][CoordinateX] = UserMatrix[0][CoordinateX];
	}

	AxillaryVector.NullMatrixFilling();

	SubstractionVector.UnitVectorFilling();

	double Coefficient = 0;

	int Iterator = 0;

	cout << "System of linear equations :\n" << endl;

	cout << UserMatrix << endl;

	while (VectorRate(SubstractionVector) > PrecisionOfResult)
	{
		Matrix TemporaryVector(1, UserMatrix.GetLengthOfMatrix() - 1);

		for (int CoordinateX = 0; CoordinateX < TemporaryVector.GetLengthOfMatrix(); CoordinateX++)
		{
			TemporaryVector[0][CoordinateX] = UserMatrix[Iterator][CoordinateX];
		}

		Coefficient = (UserMatrix[Iterator][UserMatrix.GetLengthOfMatrix() - 1] - ScalarProduct(TemporaryVector, Vector)) / (VectorRate(TemporaryVector) * VectorRate(TemporaryVector));

		TemporaryVector = TemporaryVector * Coefficient;

		AxillaryVector = Vector + TemporaryVector;

		SubstractionVector = AxillaryVector - Vector;

		Vector = AxillaryVector;

		TemporaryVector.Dispose();

		if (Iterator < TemporaryVector.GetLengthOfMatrix() - 1)
		{
			Iterator++;
		}
		else
		{
			Iterator = 0;
		}
	}

	cout << "Vector of solutions :\n" << endl;

	cout << Vector << endl;

	cout << "Last coefficient :\n" << endl;

	cout << Coefficient << endl;

	cout << endl;
}

void JakobiRotationMethod(Matrix UserMatrix, const double PrecisionOfResult)
{
	Matrix CopiedMatrix = UserMatrix,

		MatrixOfSolutions(UserMatrix.GetHeightOfMatrix()),

		BufferMatrix(UserMatrix.GetHeightOfMatrix()),

		RotationMatrix(UserMatrix.GetHeightOfMatrix()),

		ValidationMatrix(UserMatrix.GetHeightOfMatrix());

	int MaxCoordinateX = 0,

		MaxCoordinateY = 0,

		CoordinateX = 0,

		CoordinateY = 0,

		SizeOfSquareMatrix = UserMatrix.GetHeightOfMatrix();

	double MaxCoefficient = 0,

		AngleFi = 0,

		Fault = 0;

	MatrixOfSolutions.UnitMatrixFilling();

	cout << "System of linear equations :\n" << endl;

	cout << CopiedMatrix << endl;

	for (int CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
	{
		for (int CoordinateY = CoordinateX + 1; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
		{
			Fault = Fault + CopiedMatrix[CoordinateX][CoordinateY] * CopiedMatrix[CoordinateX][CoordinateY];
		}
	}

	Fault = sqrt(2 * Fault);

	while (Fault > PrecisionOfResult)
	{
		MaxCoefficient = 0;

		for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
		{
			for (CoordinateY = CoordinateX + 1; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
			{
				if (CopiedMatrix[CoordinateX][CoordinateY] > 0 && CopiedMatrix[CoordinateX][CoordinateY] > MaxCoefficient)
				{
					MaxCoefficient = CopiedMatrix[CoordinateX][CoordinateY];

					MaxCoordinateX = CoordinateX;

					MaxCoordinateY = CoordinateY;
				}
				else if (CopiedMatrix[CoordinateX][CoordinateY] < 0 && -CopiedMatrix[CoordinateX][CoordinateY] > MaxCoefficient)
				{
					MaxCoefficient = -CopiedMatrix[CoordinateX][CoordinateY];

					MaxCoordinateX = CoordinateX;

					MaxCoordinateY = CoordinateY;
				}
			}
		}

		RotationMatrix.UnitMatrixFilling();

		if (CopiedMatrix[MaxCoordinateX][MaxCoordinateX] == CopiedMatrix[MaxCoordinateY][MaxCoordinateY])
		{
			RotationMatrix[MaxCoordinateX][MaxCoordinateX] = RotationMatrix[MaxCoordinateY][MaxCoordinateY] = RotationMatrix[MaxCoordinateY][MaxCoordinateX] = sqrt(2.0) / 2.0;

			RotationMatrix[MaxCoordinateX][MaxCoordinateY] = -sqrt(2.0) / 2.0;
		}
		else
		{
			AngleFi = 0.5 * atan((2.0 * CopiedMatrix[MaxCoordinateX][MaxCoordinateY]) / (CopiedMatrix[MaxCoordinateX][MaxCoordinateX] - CopiedMatrix[MaxCoordinateY][MaxCoordinateY]));

			RotationMatrix[MaxCoordinateX][MaxCoordinateX] = RotationMatrix[MaxCoordinateY][MaxCoordinateY] = cos(AngleFi);

			RotationMatrix[MaxCoordinateX][MaxCoordinateY] = -sin(AngleFi);

			RotationMatrix[MaxCoordinateY][MaxCoordinateX] = sin(AngleFi);
		}

		BufferMatrix.NullMatrixFilling();

		for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
		{
			for (CoordinateY = 0; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
			{
				for (int AdditionalCoordinate = 0; AdditionalCoordinate < SizeOfSquareMatrix; AdditionalCoordinate++)
				{
					BufferMatrix[CoordinateX][CoordinateY] = BufferMatrix[CoordinateX][CoordinateY] + RotationMatrix[AdditionalCoordinate][CoordinateX] * CopiedMatrix[AdditionalCoordinate][CoordinateY];
				}
			}
		}

		CopiedMatrix.NullMatrixFilling();

		CopiedMatrix = BufferMatrix * RotationMatrix;

		Fault = 0.0;

		for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
		{
			for (CoordinateY = CoordinateX + 1; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
			{
				Fault = Fault + CopiedMatrix[CoordinateX][CoordinateY] * CopiedMatrix[CoordinateX][CoordinateY];
			}
		}

		Fault = sqrt(2 * Fault);

		BufferMatrix.NullMatrixFilling();

		BufferMatrix = MatrixOfSolutions * RotationMatrix;

		MatrixOfSolutions = BufferMatrix;
	}

	cout << "A solution of current system of linear equations :\n" << endl;

	for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
	{
		cout << setw(2) << "x" << CoordinateX + 1 << " = (";

		for (CoordinateY = 0; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
		{
			cout << MatrixOfSolutions[CoordinateX][CoordinateY];

			if (CoordinateY < SizeOfSquareMatrix - 1)
			{
				cout << ",";
			}

			if (CoordinateY == SizeOfSquareMatrix - 1)
			{
				cout << ")" << endl;
			}
		}
	}

	cout << endl;

	cout << "Own values of matrix :\n" << endl;

	for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
	{
		cout << setw(2) << "v" << CoordinateX + 1 << " = " << CopiedMatrix[CoordinateX][CoordinateX] << endl;
	}

	cout << endl;
}

void LinearRegression(Matrix MatrixOfCoordinates)
{
	Matrix AxillaryMatrix(MatrixOfCoordinates.GetLengthOfMatrix(), MatrixOfCoordinates.GetLengthOfMatrix() + 1);

	int CoordinateX = 0,

		CoordinateY = 0,

		AdditionalCoordinate = 0;

	AxillaryMatrix.NullMatrixFilling();

	AxillaryMatrix[0][0] = MatrixOfCoordinates.GetHeightOfMatrix();

	cout << "A solution of current linear regression :\n" << endl;

	cout << MatrixOfCoordinates << endl;

	for (CoordinateX = 1; CoordinateX < MatrixOfCoordinates.GetLengthOfMatrix(); CoordinateX++)
	{
		for (CoordinateY = 0; CoordinateY < MatrixOfCoordinates.GetHeightOfMatrix(); CoordinateY++)
		{
			AxillaryMatrix[0][CoordinateX] += MatrixOfCoordinates[CoordinateY][CoordinateX];
		}

		AxillaryMatrix[CoordinateX][0] = AxillaryMatrix[0][CoordinateX];
	}

	for (CoordinateX = 1; CoordinateX < MatrixOfCoordinates.GetLengthOfMatrix(); CoordinateX++)
	{
		for (CoordinateY = 1; CoordinateY < MatrixOfCoordinates.GetLengthOfMatrix(); CoordinateY++)
		{
			for (AdditionalCoordinate = 0; AdditionalCoordinate < MatrixOfCoordinates.GetHeightOfMatrix(); AdditionalCoordinate++)
			{
				AxillaryMatrix[CoordinateX][CoordinateY] += MatrixOfCoordinates[AdditionalCoordinate][CoordinateX] * MatrixOfCoordinates[AdditionalCoordinate][CoordinateY];
			}
		}
	}

	cout << "Matrix :\n" << endl;

	cout << AxillaryMatrix << endl;

	for (CoordinateX = 0; CoordinateX < MatrixOfCoordinates.GetHeightOfMatrix(); CoordinateX++)
	{
		AxillaryMatrix[AxillaryMatrix.GetLengthOfMatrix()][AxillaryMatrix.GetLengthOfMatrix()] += MatrixOfCoordinates[CoordinateX][0];
	}

	for (CoordinateX = 1; CoordinateX < MatrixOfCoordinates.GetLengthOfMatrix(); CoordinateX++)
	{
		for (CoordinateY = 0; CoordinateY < MatrixOfCoordinates.GetHeightOfMatrix(); CoordinateY++)
		{
			AxillaryMatrix[CoordinateX][AxillaryMatrix.GetLengthOfMatrix()] += MatrixOfCoordinates[CoordinateY][0] * MatrixOfCoordinates[CoordinateY][CoordinateX];
		}
	}

	cout << "Vector of solutions :\n" << endl;

	cout << "x = (";

	for (CoordinateX = 0; CoordinateX < AxillaryMatrix.GetHeightOfMatrix(); CoordinateX++)
	{
		cout << AxillaryMatrix[CoordinateX][AxillaryMatrix.GetLengthOfMatrix()];

		if (CoordinateX < AxillaryMatrix.GetHeightOfMatrix - 1)
		{
			cout << ",";
		}
	}

	cout << ")" << endl;
}

int ConvertToInt(int Value)
{
	unsigned char CharOne, CharTwo, CharThree, CharFour;

	CharOne = Value & 255;

	CharTwo = (Value >> 8) & 255;

	CharThree = (Value >> 16) & 255;

	CharFour = (Value >> 24) & 255;

	return((int)CharOne << 24) + ((int)CharTwo << 16) + ((int)CharThree << 8) + CharFour;
}

void MNISTDataReader(int NumberOfImageToShow, int DataOfAnImage, vector<vector<double>> &DataBase)
{
	DataBase.resize(NumberOfImageToShow, vector<double>(DataOfAnImage));

	ifstream file("C:\\t10k-images-idx3-ubyte.gz", ios::binary);

	if (file.is_open())
	{
		int Number = 0,

			NumberOfImages = 0,

			Height = 0,

			Length = 0;

		file.read((char*)&Number, sizeof(Number));

		Number = ConvertToInt(Number);

		file.read((char*)&NumberOfImages, sizeof(NumberOfImages));

		NumberOfImages = ConvertToInt(NumberOfImages);

		file.read((char*)&Height, sizeof(Height));

		Height = ConvertToInt(Height);

		file.read((char*)&Length, sizeof(Length));

		Length = ConvertToInt(Length);

		for (int CoordinateX = 0; CoordinateX < NumberOfImages; ++CoordinateX)
		{
			for (int CoordinareY = 0; CoordinareY < Height; ++CoordinareY)
			{
				for (int AdditionalCoordinate = 0; AdditionalCoordinate < Length; ++AdditionalCoordinate)
				{
					unsigned char TemporaryValue = 0;

					file.read((char*)&TemporaryValue, sizeof(TemporaryValue));

					DataBase[CoordinateX][(Height*CoordinareY) + AdditionalCoordinate] = (double)TemporaryValue;
				}
			}
		}
	}
}

void LinearClassification()
{

}

void main()
{
	vector<vector<double>> DataBase;

	MNISTDataReader(10000, 784, DataBase);

	int Counter = 1;

	for (int CoordinateX = 0; CoordinateX < 10; CoordinateX++)
	{
		for (int CoordinateY = 0; CoordinateY < 784; CoordinateY++)
		{
			cout << setw(3) << DataBase[CoordinateX][CoordinateY] << " ";

			if (Counter == 28)
			{
				cout << endl;

				Counter = 0;
			}

			Counter++;
		}
	}

	system("pause");
}