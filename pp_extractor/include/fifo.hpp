#ifndef __GENERIC_KEEP_FIFO_HPP__
#define __GENERIC_KEEP_FIFO_HPP__

template <typename T, class Pred>
class GenericKeepFifo {
public:
    GenericKeepFifo(const size_t size, Pred p) :
        size(size),
        mid(size/2),
        p(p) {
    }

    void insert(T item) {
        items.push_back(FIFOItem(item));

        if (items.size() == size) {
            // When we reach size we need to check for predicate in front
            for (size_t i = 0; i <= mid; ++i) {
                // If predicate (e.g., small PP)
                if (p(items[mid].item)) {
                    keep();
                    break;
                }
            }
        } else if (items.size() > size) {
            // When size is passed, pop front
            items.pop_front();

            // If predicate (e.g., small PP)
            if (p(items[mid].item)) {
                // Keep the information
                keep();
            }
        }
    }

    void finalize() {
        // Search for predicate (e.g., small PP) at the end
        for (size_t i = ((items.size() < size) ? 0 : (mid + 1)); i < items.size(); ++i) {
            if (p(items[i].item)) {
                keep();
                break;
            }
        }
    }

    const std::vector<T>& get_kept_items_ref() const {
        return kept_items;
    }

private:
    class FIFOItem {
    public:
        FIFOItem(T item) : item(item), kept(false) {}
        T item;
        bool kept;
    };

    void keep() {
        for (size_t i = 0; i < items.size(); ++i) {
            if (!items[i].kept) {
                kept_items.push_back(items[i].item);
                items[i].kept = true;
            }
        }
    }

protected:
    size_t size;
    size_t mid;
    /// @brief FIFO Items
    std::deque<FIFOItem> items;
    /// @brief Kept items surrounding predicate
    std::vector<T> kept_items;
    /// @brief Predicate that tells us when to keep
    Pred p;
};

#endif /* __GENERIC_KEEP_FIFO_HPP__ */
