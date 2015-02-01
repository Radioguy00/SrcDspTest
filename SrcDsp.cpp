// SrcDsp.cpp : Defines the entry point for the console application.
//

#ifdef _WIN32
#include "stdafx.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

#ifdef _WIN32
#include "..\..\SrcLibraries\SrcDsp\filters.h"
#include "..\..\SrcLibraries\SrcDsp\generators.h"
#include "..\..\SrcLibraries\SrcDsp\files.h"
#include "..\..\SrcLibraries\SrcDsp\upsampling_filters.h"
#else
#include "../../src_libraries/dsp/filters.h"
#include "../../src_libraries/dsp/generators.h"
#include "../../src_libraries/dsp/files.h"
#include "../../src_libraries/dsp/upsampling_filters.h"
#endif


bool testFilters();
bool testGenerators();
bool testFiles();
bool testUpsamplingFilters();
int common_main();

#ifdef _WIN32
int _tmain(int argc, _TCHAR* argv[])
{
	return common_main();
}
#else
int main(int argc, char ** argv)
{
	return common_main();
}

#endif

int common_main()
{
	testFiles();
	testGenerators();
	testFilters();
	testUpsamplingFilters();

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
#ifdef CPLUPLUS11
		using valueType = double;
#else
		typedef double valueType;
#endif
		GenSine<valueType> gen(1, 34.5);
		vector<valueType> out(20);

		gen.step(out);

#ifdef CPLUSPLUS11
		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
#else
		for(size_t index = 0; index < out.size(); index++)
			cout << out[index] << ' ';
		cout << "\n\n";
#endif
	}

	{
		typedef int16_t valueType ;
		GenSine<valueType> gen(1, 1000);
		vector<valueType> out(20);

		gen.step(out);

#ifdef CPLUSPLUS11
		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
#else
		for(size_t index = 0; index < out.size(); index++)
			cout << out[index] << ' ';
		cout << "\n\n";
#endif
	}


	{
#ifdef CPLUSPLUS11
		using valueType = complex<int16_t>;
#else
		typedef complex<int16_t> valueType;
#endif
		GenSine<valueType> gen(1, 1000);
		vector<valueType> out(20);

		gen.step(out);

#ifdef CPLUSPLUS11
		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
#else
		for(size_t index = 0; index < out.size(); index++)
			cout << out[index] << ' ';
		cout << "\n\n";
#endif
	}


	return false;

}

bool testUpsamplingFilters()
{
	using namespace std;

	cout << "\n **** Upsampling Filters Test" << '\n';
	const unsigned ratio = 2; // Upsampling ratio

	{
		cout << "Floating point test" << '\n';
		cout << "Creation of filter object" << '\n';
		vector<float> filterCoeff{ 1.34f, 4.56f, 0.34f, 1.08f, 5.67f, -2.54f};
		FilterUpsamplingFir<float, double, double, float, ratio> filter(filterCoeff);

		// ---- Display of filter parameters
		cout << "Upsampling Ratio: " << ratio << '\n';
		for (size_t index = 0; index< filterCoeff.size(); ++index)
			cout << filterCoeff[index] << " ";
		cout << '\n';

		// ----- Impulse response
		cout << "Filter inpulse response" << '\n';
		vector<float> input{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		vector<double> output(input.size() * ratio);
		filter.step(input, output);
		cout << "Input:\n";
		for (size_t index = 0; index< input.size(); ++index)
			cout << input[index] << " ";
		cout << '\n';
		cout << "Output:\n";
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
		cout << '\n';

		// ----- Delayed impulse response
		cout << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size() * ratio);
		filter.step(input, output);
		cout << "Input:\n";
		for (size_t index = 0; index< input.size(); ++index)
			cout << input[index] << " ";
		cout << '\n';
		cout << "Output:\n";
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
		cout << '\n';
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

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif

		cout << '\n' << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif
		cout << '\n';
	}
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
#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif
		cout << '\n';

		cout << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif
		cout << '\n';
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

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif

		cout << '\n' << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif
		cout << '\n';
	}

	return false;
}

