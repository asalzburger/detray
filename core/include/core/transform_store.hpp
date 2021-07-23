/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray
{

    using transform3 = __plugin::transform3;

    /** A static inplementation of an alignable transform store */
    template < template <typename> class vector_type = dvector>
    class static_transform_store
    {
    public:

        /** Elementwise access. Needs []oprator for storage type for now */
        inline auto operator[](const unsigned int i) { return _data[i]; }
        inline auto operator[](const unsigned int i) const { return _data[i]; }

        /** Empty context type struct */
        struct context
        {
        };

        /** Helper struct to pass range and context */
        template<typename range_iterator>
        struct contextual_range {
            using context = static_transform_store::context;
            /*context ctx;*/

            const range_iterator r;

            inline auto begin() { return r.begin(); }
            inline auto end() { return r.end(); }

            inline auto operator[](const unsigned int i) { return *(r.begin() + i); }
            inline auto operator[](const unsigned int i) const { return *(r.begin() + i); }
        };

        using storage = vector_type<transform3>;

        /** Forward iterator : Contextual STL like API
         *
         * @param ctx The context of the call (ignored)
         */
        auto begin(const context & /*ctx*/) const ->decltype(auto)
        {
            return _data.begin();
        }

        /** Forward iterator : Contextual STL like API
         *
         * @param ctx The context of the call (ignored)
         */
        auto end(const context & /*ctx*/) const ->decltype(auto)
        {
            return _data.end();
        }

        /** Access to a predefined range of elements
         *
         * @tparam start start index of rage
         * @tparam end end index of rage
         *
         * @param ctx The context of the call (ignored)
         *
         * @return range restricted iterator
         */
        const inline auto range(const size_t begin, const size_t end, const context & ctx) const
        {
            return contextual_range<decltype(range_iter(_data, dindex_range{begin, end}))>{range_iter(_data, dindex_range{begin, end})};
        }

        /** Reserve memory : Contextual STL like API
         *
         * @param ctx The context of the call (ignored)
         * @param tidx The size of the reserved container memory
         */
        void reserve(const context & /*ctx*/, size_t n_size)
        {
            _data.reserve(n_size);
        }

        /** Push back : Contextual STL like API, copy semantics
         *
         * @param ctx The context of the call (ignored)
         * @param tcf The transform to be filled
         */
        void push_back(const context & /*ctx*/, const transform3 &tf)
        {
            _data.push_back(tf);
        }

        /** Push back : Contextual STL like API, move semantics
         *
         * @param ctx The context of the call (ignored)
         * @param tcf The transform to be filled
         */
        void push_back(const context & /*ctx*/, const transform3 &&tf)
        {
            _data.push_back(std::move(tf));
        }

        /** Size : Contextual STL like API
         * @param ctx The context of the call (ignored)
         */ 
        const size_t size(const context & /*ctx*/) const
        {
            return _data.size();
        }

        /** Empty : Contextual STL like API
         * @param ctx The context of the call (ignored)
         */
        bool empty(const context & /*ctx*/)
        {
            return _data.empty();
        }

        /** Add a new bunch of (contextual) transforms - move semantics
         *
         * @param ctx The context of the call (ignored)
         * @param trfs The transform container, move semantics
         *
         * @note in general can throw an exception
         */
        void set_contextual_transforms(const context & /*ctx*/, storage &&trfs) noexcept(false)
        {
            _data = std::move(trfs);
        }

        /** Append a bunch of (contextual) transforms - move semantics
         *
         * @param ctx The context of the call (ignored)
         * @param trfs The transform container, move semantics
         *
         * @note in general can throw an exception
         */
        void append_contextual_transforms(const context & /*ctx*/, storage &&trfs) noexcept(false)
        {
            _data.reserve(_data.size() + trfs.size());
            _data.insert(_data.end(), std::make_move_iterator(trfs.begin()), std::make_move_iterator(trfs.end()));
        }

        /** Get the contextual transform
         *
         * @param ctx The context of the call (ignored)
         * @param tidx The transform index
         */
        const transform3 &contextual_transform(const context & /*ctx*/, dindex tidx) const
        {
            return _data[tidx];
        }

    private:
        /** Common to surfaces & portals: transform store */
        vector_type<transform3> _data;
    };

} // namespace detray
