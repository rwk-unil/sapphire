#ifndef __GENERIC_KEEP_FIFO_HPP__
#define __GENERIC_KEEP_FIFO_HPP__

template <typename T, class Pred>
class GenericKeepFifo {
public:
    GenericKeepFifo(const size_t size, Pred p) :
        size(size),
        mid(size/2),
        p(p) {
        if (!(size & 1)) {
            std::cerr << "FIFO size should be odd ! Adjusting size to " << ++this->size << std::endl;
            this->mid = this->size/2;
        }
    }

    void insert(T item) {
        // If the FIFO is empty, fill with "dummy items" (to simplify logic)
        if (items.size() == 0) {
            for (size_t i = 0; i < size; ++i) {
                items.push_back(FIFOItem(item));
                // Setting this true will make sure "keep" doesn't save them
                items.back().kept = true;
            }
        }

        // Push item at the end
        items.push_back(FIFOItem(item));

        // Pop the oldest item
        items.pop_front();

        // If predicate (e.g., small PP)
        if (p(items[mid].item)) {
            // Keep the information
            keep();
        }
    }

    void finalize() {
        // Search for predicate (e.g., small PP) at the end
        for (size_t i = mid + 1; i < items.size(); ++i) {
            if (p(items[i].item)) {
                keep_end(i-mid);
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

    inline void keep_end(size_t start) {
        for (size_t i = start; i < items.size(); ++i) {
            if (!items[i].kept) {
                kept_items.push_back(items[i].item);
                items[i].kept = true;
            }
        }
    }

    inline void keep() {
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
