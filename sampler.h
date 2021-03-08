template <typename T>
class sampler
{
    std::vector<T> keys;
    std::discrete_distribution<T> distr;

public:
    sampler(const std::vector<T>& keys, const std::vector<float>& prob) :
        keys(keys), distr(prob.begin(), prob.end()) { }

    T operator()()
    {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        return keys[distr(gen)];
    }
};