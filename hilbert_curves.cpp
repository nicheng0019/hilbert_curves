#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>

std::string _binary_repr(int num, int width)
{
	std::string str;
	while (1)
	{
		int b = num & 0x1;
		str = std::to_string(b) + str;

		num = num >> 1;
		if (num == 0)
		{
			break;
		}
	}

	while (str.length() < width)
	{
		str = "0" + str;
	}

	return str;
}

class HilbertCurve
{
public:
	HilbertCurve(int p, int n)
	{
		m_p = p;
		m_n = n;
		m_min_h = 0;
		m_max_h = pow(2, m_p * m_n) - 1;
		m_min_x = 0;
		m_max_x = pow(2, m_p) - 1;
	};



	std::vector<int> _hilbert_integer_to_transpose(int h)
	{
		std::vector<int> x;
		std::string str = _binary_repr(h, m_p * m_n);
		for (int i = 0; i < m_n; i++)
		{
			std::string substr;
			for (int n = i; n < str.length(); n += m_n)
			{
				substr += str.at(n);
			}
			int v = std::stoi(substr, nullptr, 2);
			x.push_back(v);
		}

		return x;
	}

	int _transpose_to_hilbert_integer(std::vector<int> x)
	{
		std::vector<std::string> x_bit_str;
		for (int i = 0; i < m_n; i++)
		{
			std::string str = _binary_repr(x[i], m_p);
			x_bit_str.push_back(str);
		}

		std::string totalstr;
		for (int i = 0; i < m_p; i++)
		{
			for (auto str : x_bit_str)
			{
				totalstr += str[i];
			}
		}
	
		return std::stoi(totalstr, nullptr, 2);
	}

	std::vector<int> point_from_distance(int distance)
	{
		std::vector<int> x = _hilbert_integer_to_transpose(distance);
		int z = 2 << (m_p - 1);

		int t = x[m_n - 1] >> 1;
		for (int i = m_n - 1; i > 0; i--)
		{
			x[i] ^= x[i - 1];
		}

		x[0] ^= t;

		int q = 2;
		while (q != z)
		{
			int p = q - 1;
			for (int i = m_n - 1; i > -1; i--)
			{
				if (x[i] & q)
				{
					x[0] ^= p;
				}
				else
				{
					t = (x[0] ^ x[i]) & p;
					x[0] ^= t;
					x[i] ^= t;
				}
			}

			q = q << 1;
		}

		return x;
	}

	std::vector<std::vector<int>> points_from_distances(std::vector<int> distances)
	{
		std::vector<std::vector<int>> points;
		for (auto distance : distances)
		{
			std::vector<int> x = point_from_distance(distance);
			points.push_back(x);
		}

		return points;
	}

	int distance_from_point(std::vector<int> point)
	{
		int m = 1 << (m_p - 1);
		int q = m;

		while (q > 1)
		{
			int p = q - 1;
			for (int i = 0; i < m_n; i++)
			{
				if (point[i] & q)
					point[0] ^= p;
				else
				{
					int t = (point[0] ^ point[i]) & p;
					point[0] ^= t;
					point[i] ^= t;
				}
			}

			q = q >> 1;
		}

		for (int i = 1; i < m_n; i++)
		{
			point[i] ^= point[i - 1];
		}

		int t = 0;
		q = m;

		while (q > 1)
		{
			if (point[m_n - 1] & q)
			{
				t ^= q - 1;
			}
			q = q >> 1;
		}

		for (int i = 0; i < m_n; i++)
		{
			point[i] ^= t;
		}

		return _transpose_to_hilbert_integer(point);
	}

	std::vector<int> distances_from_points(std::vector<std::vector<int>> points)
	{
		std::vector<int> distances;
		for (auto& point : points)
		{
			int distance = distance_from_point(point);
			distances.push_back(distance);
		}

		return distances;
	}

	int p()
	{
		return m_p;
	}

	int n()
	{
		return m_n;
	}

private:
	int m_p;
	int m_n;
	int m_min_h;
	int m_max_h;
	int m_min_x;
	int m_max_x;
};


void pcnormalize(std::vector<std::vector<double>> xyz, std::vector<std::vector<double>>& xyz_normelised, std::vector<std::vector<double>>& Limits)
{
	int count = xyz.size();
	int n = xyz[0].size();

	for (int i = 0; i < n; i++)
	{
		std::vector <double> ndata;
		for (int j = 0; j < count; j++)
		{
			ndata.push_back(xyz[j][i]);
		}
		std::pair< std::vector <double>::iterator, std::vector <double>::iterator >Position = std::minmax_element(ndata.begin(), ndata.end());
		double minv = *(Position.first);
		double maxv = *(Position.second);

		Limits.push_back({ maxv , minv });
	}

	for (int i = 0; i < count; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (Limits[j][0] - Limits[j][1] < 0.0000001)
			{
				continue;
			}
			xyz_normelised[i][j] = (xyz[i][j] - Limits[j][1]) / (Limits[j][0] - Limits[j][1]);
		}
	}
}

std::vector<int> Coordinate2Distance(std::vector<std::vector<double>> xyz_location, HilbertCurve hilbertcurve)
{
	int count = xyz_location.size();
	int n = xyz_location[0].size();

	std::vector<int> distances;
	for (int i = 0; i < count; i++)
	{
		std::vector<int> coordinate;
		for (int j = 0; j < n; j++)
		{
			int x = floor(xyz_location[i][j] * pow(2, hilbertcurve.p() - 1));
			coordinate.push_back(x);
		}
	
		int distance = hilbertcurve.distance_from_point(coordinate);

		distances.push_back(distance);
	}

	return distances;
}

template<class BidiIter >
BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random) {
	size_t left = std::distance(begin, end);
	srand(time(NULL));
	while (num_random--) {
		BidiIter r = begin;
		std::advance(r, rand() % left);
		std::swap(*begin, *r);
		++begin;
		--left;
	}
	return begin;
}

std::vector<std::vector<double>> downsample(std::vector<std::vector<double>> xyz_location, 
									std::vector<std::vector<double>> Limits, 
									std::vector<int> distances,
									double b, int subsample_num, int P)
{
	int count = xyz_location.size();
	int n = xyz_location[0].size();
	int num = floor(pow(2, (n * P)) / b) + 1;

	std::vector<std::pair<int, std::vector<int>>> bin;
	for (int i = 0; i < num; i++)
	{
		bin.push_back(std::make_pair<int, std::vector<int>>(0, {}));
	}

	for (int i = 0; i < count; i++)
	{
		int location = floor(distances[i] / b + 1);
		bin[location].first += 1;
		bin[location].second.push_back(i);
	}
	
	int non_zero = 0;
	for (int i = 0; i < num; i++)
	{
		if (bin[i].first)
		{
			non_zero++;
		}
	}

	int subsample_num_bin = ceil(subsample_num * 1.0 / non_zero);

	for (int i = 0; i < num; i++)
	{
		if (bin[i].first > subsample_num_bin)
		{
			std::vector<int> y;
			std::vector<int>::iterator end = random_unique(bin[i].second.begin(), bin[i].second.end(), subsample_num_bin);
			y.insert(y.begin(), bin[i].second.begin(), end);
			bin[i].second = y;
		}
	}

	std::vector<int> selection;
	for (int i = 0; i < num; i++)
	{
		selection.insert(selection.end(), bin[i].second.begin(), bin[i].second.end());
	}

	int ss = selection.size();

	std::vector<std::vector<double>> data;
	for (int i = 0; i < ss; i++)
	{
		data.push_back(xyz_location[selection[i]]);
		for (int j = 0; j < n; j++)
		{
			data[i][j] = data[i][j] * (Limits[j][0] - Limits[j][1]) + Limits[j][1];
		}
	}

	return data;
}

std::vector<std::vector<double>> HilbertCurve_Subsample(int P, double bin_size, int q, std::vector<std::vector<double>> xyz)
{
	int n = xyz[0].size();
	HilbertCurve hilbertCurve(P, n);

	std::vector<std::vector<double>> xyz_normalized = xyz;
	std::vector<std::vector<double>> Limits;
	pcnormalize(xyz, xyz_normalized, Limits);

	std::vector<int> distances = Coordinate2Distance(xyz_normalized, hilbertCurve);

	std::vector<std::vector<double>> xyz_subsampled = downsample(xyz_normalized, Limits, distances, bin_size, q, P);
	return xyz_subsampled;
}


int main2()
{
	FILE* f = fopen("D:/Program/Project/GNN/xyz.bin", "rb");
	double data[20 * 3];
	fread(data, sizeof(double), 20 * 3, f);
	fclose(f);

	std::vector<std::vector<double>> xyz;
	for (int i = 0; i < 20; i++)
	{
		std::vector<double> temp_data;
		for (int j = 0; j < 3; j++)
		{
			temp_data.push_back(data[i * 3 + j]);
		}
		xyz.push_back(temp_data);
	}

	int P = 2;
	double  bin_size = 0.5;
	int q = 5;
	std::vector<std::vector<double>> xyz_subsampled = HilbertCurve_Subsample(P, bin_size, q, xyz);

	return 1;
}