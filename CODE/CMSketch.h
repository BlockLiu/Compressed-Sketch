#include "BOBHash32.h"
#include "Basic.h"

using std::min;
using std::swap;

template<uint8_t key_len, Compress_Method compress_method, int d = 4>
struct CMSketch {
	int mem_in_bytes;
	int w = 1, k = 0, r = 1;
	int *cm_sketch[32][d];
	BOBHash32 *hash[d];

	public:
	string name;

	CMSketch(int mem_in_bytes_) : mem_in_bytes(mem_in_bytes_) {
		for (; w <= mem_in_bytes / 4 / d / 2; w *= 2, k += 1);

		for (int i = 0; i < d; i++) {
			if (compress_method != Hierarchical) {
				cm_sketch[0][i] = new int[w];
				memset(cm_sketch[0][i], 0, sizeof(int) * w);
			} else {
				for (int l = 1; l < k; l++) {
					cm_sketch[l][i] = new int[w >> l];
					memset(cm_sketch[l][i], 0, sizeof(int) * (w >> l));
				}
			}
		}

		random_device rd;
		for (int i = 0; i < d; i++) {
			hash[i] = new BOBHash32(uint32_t(rd() % MAX_PRIME32));
		}

		stringstream name_buf;
		name_buf << "CMSketch@" << mem_in_bytes;
		name = name_buf.str();
	}

	void insert(uint8_t *key) {
		for (int i = 0; i < d; i++) {
			if (compress_method != Hierarchical) {
				int idx = hash[i]->run((char *)key, key_len) >> (32 - k);
				cm_sketch[0][i][idx] += 1;
			} else {
				for (int l = 1; l < k; l++) {
					int idx = hash[i]->run((char *)key, key_len) >> (32 - k + l);
					cm_sketch[l][i][idx] += 1;
				}
			}
		}
	}

	int query(uint8_t *key) {
		int tmin = 1 << 30, ans = tmin;
		for (int i = 0; i < d; i++) {
			if (compress_method != Hierarchical) {
				int idx = hash[i]->run((char *)key, key_len) >> (32 - k);
				int val = cm_sketch[0][i][idx];
				ans = std::min(val, ans);
			} else {
				int idx = hash[i]->run((char *)key, key_len) >> (32 - k + r);
				int val = cm_sketch[r][i][idx];
				ans = std::min(val, ans);
			}
		}
		return ans;
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
								tmax += cm_sketch[0][i][idx + l];
								break;

							case Max_Selectionp:
								tmax = max(tmax, cm_sketch[0][i][idx + l]);
								break;
						}

					cm_sketch[0][i][j] = tmax;
				}
		} else {
			r += rate;
		}
	}

	~CMSketch() {
		for (int i = 0; i < d; i++) {
			delete hash[i];
			if (compress_method != Hierarchical) {
				delete cm_sketch[0][i];
			} else {
				for (int l = 1; l < k; l++)
					delete cm_sketch[l][i];
			}
		}
	}
};
