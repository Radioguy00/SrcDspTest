// SrcDsp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
//#include "D:\SoftApp\SrcLibraries\SrcDsp\filters.h"
#include "..\..\SrcLibraries\SrcDsp\filters.h"
#include "..\..\SrcLibraries\SrcDsp\generators.h"
#include "..\..\SrcLibraries\SrcDsp\files.h"

bool testFilters();
bool testGenerators();
bool testFiles();

int _tmain(int argc, _TCHAR* argv[])
{
	testFiles();

	return 0;
}

bool testFiles()
{
	using namespace dsptl;
	using namespace std;

	{
		std::vector<int16_t> in{ 12, 34, 56, 78, 1, 43 };
		ofstream os("file_int16.txt");
		saveAsciiSamples(in, os);
		os.close();
	}

	{
		std::vector<float> in{ 12.34f, 34.5f, 56.67f, 78.f, 1.45f, 43.f };
		ofstream os("file_float.txt");
		saveAsciiSamples(in, os);
		os.close();
	}

	{
		std::vector<int16_t> in{ 12, 34, 56, 78, 1, 43 };
		ofstream os("file_int16.bin");
		saveBinarySamples(in, os);
		os.close();
	}


	{
		std::vector<complex<int16_t> > in{ { 12, 2 }, { 34, 56 }, { 78, 1 }, { 43, -14 } };
		ofstream os("file_complexint16.txt");
		saveAsciiSamples(in, os);
		os.close();
	}

	return false;
	

}

bool testGenerators()
{
	using namespace dsptl;
	using namespace std;

	{
		using valueType = double;
		GenSine<valueType> gen(1, 34.5);
		vector<valueType> out(20);

		gen.step(out);

		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
	}

	{
		using valueType = int16_t;
		GenSine<valueType> gen(1, 1000);
		vector<valueType> out(20);

		gen.step(out);

		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
	}

	{
		using valueType = complex<int16_t>;
		GenSine<valueType> gen(1, 1000);
		vector<valueType> out(20);

		gen.step(out);

		for (const auto& elt : out)
			cout << elt.real() << ' ' << elt.imag() << '\n';
		cout << "\n\n";
	}

	return false;

}

bool testFilters()
{
	using namespace std;

	cout << "Program is starting" << '\n';
	cout << "Floating point test" << '\n';

	{
		cout << "Creation of filter object" << '\n';
		vector<float> filterCoeff{ 1.34f, 4.56f, 0.34f };
		FilterFir<float, double, double, float> filter(filterCoeff);

		cout << "Filter inpulse response" << '\n';
		vector<float> input{ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		vector<double> output(input.size());
		filter.step(input, output);

		for (const auto& elt : output)
			cout << elt << " ";

		cout << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);

		for (const auto& elt : output)
			cout << elt << " ";
	}

	cout << '\n' << "Integer test" << '\n';
	{
		cout << "Creation of filter object" << '\n';
		vector<int16_t> filterCoeff{ 234, -23, 12 };
		FilterFir<int16_t, int32_t, int32_t, int16_t> filter(filterCoeff);

		cout << "Filter inpulse response" << '\n';
		vector<int16_t> input{ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		vector<int32_t> output(input.size());
		filter.step(input, output);

		for (const auto& elt : output)
			cout << elt << " ";

		cout << '\n' << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);

		for (const auto& elt : output)
			cout << elt << " ";
	}

	return false;
}