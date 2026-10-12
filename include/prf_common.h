#pragma once

#include <volePSI/Paxos.h>

// Aligned.h requires the span/u64 declarations supplied by Paxos.h's includes.
#include <cryptoTools/Common/Aligned.h>
#include <limits>
#include <memory>
#include <stdexcept>

namespace fpsi {
    using oc::block;
    using oc::Matrix;
    using oc::MatrixView;
    using oc::PRNG;
    using oc::span;
    using oc::u64;
    using oc::u8;
    using volePSI::Baxos;
    using volePSI::PaxosParam;
    using Socket = coproto::Socket;
    using Proto = coproto::task<void>;

    namespace detail {
        // Check before multiplying, allocating or constructing a matrix/span.
        inline u64 checkedBytes(u64 rows, u64 width)
        {
            const auto limit = static_cast<u64>(std::numeric_limits<std::ptrdiff_t>::max());
            if (width && rows > limit / width)
                throw std::length_error("PRF buffer size exceeds addressable memory");
            return rows * width;
        }

        inline u64 chunkElements(u64 bytes, u64 elementBytes)
        {
            if (!bytes || bytes >= std::numeric_limits<oc::u32>::max() || bytes % elementBytes)
                throw std::invalid_argument("PRF chunk must be nonzero, element-aligned and smaller than UINT32_MAX bytes");
            return bytes / elementBytes;
        }

        // coproto may finish a moved send before its bytes leave the socket. Each
        // queued slice owns the allocation, including when a coroutine exits by exception.
        template <class T>
        struct SharedSlice {
            using value_type = T;
            using size_type = std::size_t;
            std::shared_ptr<oc::AlignedUnVector<T>> owner;
            u64 offset = 0, count = 0;
            T *data()
            {
                return owner->data() + offset;
            }
            const T *data() const
            {
                return owner->data() + offset;
            }
            std::size_t size() const
            {
                return count;
            }
        };
    } // namespace detail
} // namespace fpsi
