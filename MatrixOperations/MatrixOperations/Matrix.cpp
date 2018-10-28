#define _USE_MATH_DEFINES
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

	double* operator [] (int Index)
	{
		return SquareMatrix[Index];
	}

	Matrix operator ~ ()
	{
		double Buffer = 0;

		int Coefficient = 1;

		if (HeightOfMatrix == LengthOfMatrix)
		{
			Coefficient = 0;
		}

		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix - Coefficient; CoordinateY++)
			{
				Buffer = SquareMatrix[CoordinateX][CoordinateY];

				SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateY][CoordinateX];

				SquareMatrix[CoordinateY][CoordinateX] = Buffer;
			}
		}

		return *this;
	}

	Matrix operator * (double Multiplicator)
	{
		for (int CoordinateX = 0; CoordinateX < HeightOfMatrix; CoordinateX++)
		{
			for (int CoordinateY = 0; CoordinateY < LengthOfMatrix; CoordinateY++)
			{
				SquareMatrix[CoordinateX][CoordinateY] = SquareMatrix[CoordinateX][CoordinateY] * Multiplicator;
			}
		}

		return *this;
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
					MultiplicationMatrix[CoordinateX][CoordinateY] = 0;

					for (int AdditionalCoordinate = 0; AdditionalCoordinate < SecondMatrix.HeightOfMatrix; AdditionalCoordinate++)
					{
						MultiplicationMatrix.SquareMatrix[CoordinateX][CoordinateY] += SquareMatrix[CoordinateX][AdditionalCoordinate] * SecondMatrix.SquareMatrix[AdditionalCoordinate][CoordinateY];
					}
				}
			}

			return MultiplicationMatrix;
		}
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
		cout << "System of linear equations :\n" << endl;

		if (UserMatrix.HeightOfMatrix != UserMatrix.LengthOfMatrix)
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
		else
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

		return OutputStream;
	}

	friend istream& operator >> (istream& InputStream, const Matrix& UserMatrix)
	{
		if (UserMatrix.HeightOfMatrix != UserMatrix.LengthOfMatrix)
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
		else
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

double ScalarProduct(Matrix UserMatrix)
{
	double Scalar = 0;

	for (int CoordinateX = 0; CoordinateX < UserMatrix.GetLengthOfMatrix(); CoordinateX++)
	{
		Scalar += UserMatrix[0][CoordinateX] * UserMatrix[0][CoordinateX];
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
	Matrix Vector1(1, UserMatrix.GetLengthOfMatrix());

	Matrix Vector2(1, UserMatrix.GetLengthOfMatrix());

	Matrix Vector3(1, UserMatrix.GetLengthOfMatrix());

	int Counter = 0;

	for (int CoordinateX = 0; CoordinateX < UserMatrix.GetLengthOfMatrix(); CoordinateX++)
	{
		Vector1[0][CoordinateX] = UserMatrix[0][CoordinateX];
	}

	Vector2.NullMatrixFilling();

	Vector3.UnitVectorFilling();

	while (sqrt(ScalarProduct(Vector3)) > PrecisionOfResult)
	{
		Matrix Vector4(1,UserMatrix.GetLengthOfMatrix());

		for (int CoordinateX = 0; CoordinateX < UserMatrix.GetLengthOfMatrix(); CoordinateX++)
		{
			UserMatrix[Counter][CoordinateX] = Vector4[0][CoordinateX];
		}

		double Value = (UserMatrix[0][Counter] - ScalarProduct(Vector1)) / (sqrt(ScalarProduct(Vector4))*sqrt(ScalarProduct(Vector4)));

		Vector4 = Vector4 * Value;

		Vector2 = Vector1 + Vector4;

		Vector3 = Vector2 - Vector1;

		Vector1 = Vector2;

		if (Counter < UserMatrix.GetLengthOfMatrix())
		{
			Counter++;
		}
		else
		{
			Counter = 0;
		}
	}

	cout << Vector1 << endl;
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

		for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
		{
			for (CoordinateY = 0; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
			{
				for (int AdditionalCoordinate = 0; AdditionalCoordinate < SizeOfSquareMatrix; AdditionalCoordinate++)
				{
					CopiedMatrix[CoordinateX][CoordinateY] = CopiedMatrix[CoordinateX][CoordinateY] + BufferMatrix[CoordinateX][AdditionalCoordinate] * RotationMatrix[AdditionalCoordinate][CoordinateY];
				}
			}
		}

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

		for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
		{
			for (CoordinateY = 0; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
			{
				for (int AdditionalCoordinate = 0; AdditionalCoordinate < SizeOfSquareMatrix; AdditionalCoordinate++)
				{
					BufferMatrix[CoordinateX][CoordinateY] = BufferMatrix[CoordinateX][CoordinateY] + MatrixOfSolutions[CoordinateX][AdditionalCoordinate] * RotationMatrix[AdditionalCoordinate][CoordinateY];
				}
			}
		}

		MatrixOfSolutions = BufferMatrix;
	}

	cout << "A solution of current system of linear equations :\n" << endl;

	for (CoordinateX = 0; CoordinateX < SizeOfSquareMatrix; CoordinateX++)
	{
		cout << setw(2) << "x" << CoordinateX + 1 << " = (";

		for (CoordinateY = 0; CoordinateY < SizeOfSquareMatrix; CoordinateY++)
		{
			cout << MatrixOfSolutions[CoordinateY][CoordinateX];

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

void main()
{
	Matrix UserMatrix(1, 4);

	UserMatrix.AsimmetricalMatrixFilling();

	KaczmarzMethod(UserMatrix, 0.0000001);

	system("pause");
}