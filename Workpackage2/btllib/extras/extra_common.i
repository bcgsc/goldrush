%ignore operator<<;

%ignore btllib::DataStream;
%ignore btllib::DataStream::operator FILE*() const;
%ignore btllib::DataSource::operator FILE*() const;
%ignore btllib::DataSink::operator FILE*() const;

%ignore btllib::CString;
%ignore btllib::CString::operator=;

%ignore btllib::OrderQueue;
%ignore btllib::OrderQueue::Block;
%ignore btllib::OrderQueue::Slot;
%ignore btllib::OrderQueue::Block::operator=;
%ignore btllib::OrderQueue::Slot::operator=;

%ignore btllib::SeqReaderFastaModule;
%ignore btllib::SeqReaderFastqModule;
%ignore btllib::SeqReaderGfa2Module;
%ignore btllib::SeqReaderSamModule;
%ignore btllib::SeqReaderMultilineFastaModule;
%ignore btllib::SeqReaderMultilineFastqModule;

%rename (SeqReaderRecord) btllib::SeqReader::Record;
%rename (SeqReaderFlag) btllib::SeqReader::Flag;
%rename (IndexlrRecord) btllib::Indexlr::Record;
%rename (IndexlrFlag) btllib::Indexlr::Flag;

%ignore btllib::SeqReader::read_block;

%ignore btllib::IORedirection;
%ignore btllib::IORedirection::in;
%ignore btllib::IORedirection::out;
%ignore btllib::IORedirection::err;

%ignore btllib::ProcessPipeline;
%ignore btllib::ProcessPipeline::in;
%ignore btllib::ProcessPipeline::out;

%ignore btllib::ProcessPipelineInternal;

%ignore btllib::BLOOM_FILTER_MAGIC_HEADER;
%ignore btllib::COUNTING_BLOOM_FILTER_MAGIC_HEADER;

%ignore btllib::SeqReader::RecordIterator::operator++;
%ignore btllib::SeqReader::RecordIterator::operator!=;
%ignore btllib::SeqReader::RecordIterator::operator*;

%ignore btllib::NtHash::NtHash(const char*, size_t, unsigned, unsigned, size_t pos = 0);

%ignore btllib::BloomFilterInitializer;
%ignore btllib::BloomFilterInitializer::operator=;

%template(VectorString) std::vector<std::string>;
%template(VectorInt) std::vector<int>;
%template(VectorUnsigned) std::vector<unsigned>;
%template(VectorDouble) std::vector<double>;
%template(VectorUint64t) std::vector<uint64_t>;

namespace btllib {
using SpacedSeed = std::vector<unsigned>;
}
%template(VectorSpacedSeed) std::vector<btllib::SpacedSeed>;

%{
  #include <map>
  #include <mutex>
  static long nthash_last_id = 0;
  static std::mutex nthash_mutex;
  static std::map<long, std::string> nthash_strings;
  static std::map<void*, long> nthash_ids;
%}

%ignore btllib::NtHash::NtHash(const std::string&, unsigned, unsigned, size_t pos = 0);
%ignore btllib::NtHash::NtHash(const char*, size_t, unsigned, unsigned, size_t pos = 0);

%ignore btllib::SeedNtHash::SeedNtHash(const char*, size_t, const std::vector<SpacedSeed>&, unsigned, unsigned, size_t pos = 0);
%ignore btllib::SeedNtHash::SeedNtHash(const std::string&, const std::vector<SpacedSeed>&, unsigned, unsigned, size_t pos = 0);
%ignore btllib::SeedNtHash::SeedNtHash(const char*, size_t, const std::vector<std::string>&, unsigned, unsigned, size_t pos = 0);
%ignore btllib::SeedNtHash::SeedNtHash(const std::string&, const std::vector<std::string>&, unsigned, unsigned, size_t pos = 0);

%extend btllib::NtHash {
  NtHash(std::string seq, unsigned hash_num, unsigned k, size_t pos = 0)
  {
    std::unique_lock<std::mutex> lock(nthash_mutex);
    nthash_strings[++nthash_last_id] = std::move(seq);
    auto *nthash = new btllib::NtHash(nthash_strings[nthash_last_id], hash_num, k);
    nthash_ids[(void*)nthash] = nthash_last_id;
    return nthash;
  }

  ~NtHash() {
    std::unique_lock<std::mutex> lock(nthash_mutex);
    nthash_strings.erase(nthash_ids[(void*)self]);
    nthash_ids.erase((void*)self);
  }
}

%extend btllib::SeedNtHash {
  SeedNtHash(std::string seq, const std::vector<SpacedSeed>& seeds, unsigned hash_num_per_seed, unsigned k, size_t pos = 0)
  {
    std::unique_lock<std::mutex> lock(nthash_mutex);
    nthash_strings[++nthash_last_id] = std::move(seq);
    auto *nthash = new btllib::SeedNtHash(nthash_strings[nthash_last_id], seeds, hash_num_per_seed, k);
    nthash_ids[(void*)nthash] = nthash_last_id;
    return nthash;
  }

  ~SeedNtHash() {
    std::unique_lock<std::mutex> lock(nthash_mutex);
    nthash_strings.erase(nthash_ids[(void*)self]);
    nthash_ids.erase((void*)self);
  }
}