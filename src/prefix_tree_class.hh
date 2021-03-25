#ifndef N_FOLD_PREFIX_TREE_CLASS_HH
#define N_FOLD_PREFIX_TREE_CLASS_HH

#include <iostream>

#include "utils.hh"

template <typename U>
class prefix_tree {
  private:
    struct node {
        U val;
        size_t parentIdx;
        int childrenCnt;

        node(U _val, size_t pIdx) : val{_val}, parentIdx{pIdx}, childrenCnt{U(0)} {}
    };

  public:
    const static size_t NO_PARENT = std::numeric_limits<size_t>::max();

    std::vector<node> tree;
    std::vector<size_t> freeIndices;

    int add(U value, size_t parInd) {
        assert(tree.size() > parInd || parInd == NO_PARENT);

        if(parInd != NO_PARENT) {
            tree[parInd].childrenCnt++;
        }

        if(freeIndices.size()) {
            int idx = freeIndices.back();
            freeIndices.pop_back();
            tree[idx] = node(value, parInd);
            return idx;
        } else {
            tree.emplace_back(value, parInd);
            return tree.size() - 1;
        }
    }

    void remove(size_t index) {
        assert(tree.size() > index);
        freeIndices.push_back(index);
        if(tree[index].parentIdx != NO_PARENT) {
            int parIdx = tree[index].parentIdx;
            if(--tree[parIdx].childrenCnt == 0) {
                remove(parIdx);
            }
        }
    }

    std::vector<U> constructPath(size_t index) {
        assert(index < tree.size());

        std::vector<U> res;
        while(tree[index].parentIdx != NO_PARENT) {
            res.push_back(tree[index].val);
            index = tree[index].parentIdx;
        }
        std::reverse(res.begin(), res.end());

        return res;
    }

    void clear() {
        tree.clear();
        freeIndices.clear();
    }
};

#endif //N_FOLD_PREFIX_TREE_CLASS_HH