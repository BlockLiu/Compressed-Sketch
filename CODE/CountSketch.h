#include "BOBHash32.h"
#include "Basic.h"

using std::min;
using std::swap;

template<uint8_t key_len, Compress_Method compress_method, int d = 4>
struct CountSketch {
	int mem_in_bytes;
	int w = 1, k = 0, r = 1;
	int *count_sketch[32][d];
	BOBHash32 *hash[d], *hash_polar[d];

	public:
	string name;

	CountSketch(int mem_in_bytes_) : mem_in_bytes(mem_in_bytes_) {
		for (; w <= mem_in_bytes / 4 / d / 2; w *= 2, k += 1);

		for (int i = 0; i < d; i++) {
			if (compress_method != Hierarchical) {
				count_sketch[0][i] = new int[w];
				memset(count_sketch[0][i], 0, sizeof(int) * w);
			} else {
				for (int l = 1; l < k; l++) {
					count_sketch[l][i] = new int[w >> l];
					memset(count_sketch[l][i], 0, sizeof(int) * (w >> l));
				}
			}
		}

		random_device rd;
		for (int i = 0; i < d; i++) {
			hash[i] = new BOBHash32(uint32_t(rd() % MAX_PRIME32));
			hash_polar[i] = new BOBHash32(uint32_t(rd() % MAX_PRIME32));
		}

		stringstream name_buf;
		name_buf << "CountSketch@" << mem_in_bytes;
		name = name_buf.str();
	}

	void insert(uint8_t *key) {
		for (int i = 0; i < d; i++) {
			if (compress_method != Hierarchical) {
				int idx = hash[i]->run((char *)key, key_len) >> (32 - k);
				int polar = hash_polar[i]->run((char *)key, key_len) % 2;

				count_sketch[0][i][idx] += polar ? 1 : -1;
			} else {
				for (int l = 1; l < k; l++) {
					int idx = hash[i]->run((char *)key, key_len) >> (32 - k + l);
					int polar = hash_polar[i]->run((char *)key, key_len) % 2;

					count_sketch[l][i][idx] += polar ? 1 : -1;
				}
			}
		}
	}

	int query(uint8_t *key) {
		int tmin, ans[d];
		for (int i = 0; i < d; i++) {
			if (compress_method != Hierarchical) {
				int idx = hash[i]->run((char *)key, key_len) >> (32 - k);
				int polar = hash_polar[i]->run((char *)key, key_len) % 2;

				int val = count_sketch[0][i][idx];
				ans[i] = polar ? val : -val;
			} else {
				int idx = hash[i]->run((char *)key, key_len) >> (32 - k + r);
				int polar = hash_polar[i]->run((char *)key, key_len) % 2;

				int val = count_sketch[r][i][idx];
				ans[i] = polar ? val : -val;
			}
		}

		sort(ans, ans + d);

		if (d % 2 == 0)
		    tmin = (ans[d / 2] + ans[d / 2 - 1]) / 2;
		else
		    tmin = ans[d / 2];
		tmin = (tmin <= 1) ? 1 : tmin;

		return tmin;
	}

	int memory_use() {
		if (compress_method != Hierarchical)
			return k + 4;
		else
			return k - r + 4;
	}

	void compress(int rate) {
		if (compress_method != Hierarchical) {
			w >>= rate, k -= rate;

			for (int i = 0; i < d; i++)
				for (int j = 0; j < w; j++) {
					int idx = j << rate, tmax = 0;

					for (int l = 0; l < 1 << rate; l++)
						switch (compress_method) {
							case Sum_Selection:
								tmax += count_sketch[0][i][idx + l];
								break;

							case Max_Selectionp:
								tmax = max(tmax, count_sketch[0][i][idx + l]);
								break;
						}

					count_sketch[0][i][j] = tmax;
				}
		} else {
			r += rate;
		}
	}

	~CountSketch() {
		for (int i = 0; i < d; i++) {
			delete hash[i];
			delete hash_polar[i];
			if (compress_method != Hierarchical) {
				delete count_sketch[0][i];
			} else {
				for (int l = 1; l < k; l++)
					delete count_sketch[l][i];
			}
		}
	}
};
