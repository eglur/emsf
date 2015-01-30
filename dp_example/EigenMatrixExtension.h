Matrix(std::string filename)
{
    loadFromFile(filename);
}


bool saveToFile(const std::string &filename, const int lim_rows = -1, const int lim_cols = -1, const int precision = 20) const
{
    std::ofstream ofs(filename.c_str());
    if (!ofs.good()) return false;

    size_t pre = ofs.precision();
    ofs.precision(precision);

    int lr = lim_rows > 0 ? lim_rows : this->rows();
    int lc = lim_cols > 0 ? lim_cols : this->cols();

    for (int i = 0; i < lr; ++i)
    {
        ofs << this->operator()(i,0);
        for (int j = 1; j < lc; ++j) ofs  << " " <<  this->operator()(i,j);
        ofs << std::endl;
    }
    ofs.close();
    ofs.precision(pre);
    return true;
}


bool loadFromFile(const std::string &filename, int nrows = -1, int ncols = -1 )
{
    if (nrows < 0 || ncols < 0)
    {

        std::ifstream ifs(filename.c_str());
        if (!ifs.good()) return false;

        std::string stmp;
        getline(ifs,stmp);
        ncols = 1;
        bool already_counted = true; // to skip initial spaces
        for (size_t i=0; i < stmp.length(); ++i)
        {
            if (stmp[i] == ' ')
            {
                if (!already_counted)
                {
                    ++ncols;
                    already_counted = true;
                }
            }
            else already_counted = false;
        }

        nrows = 0;
        while (!ifs.eof())
        {
            getline(ifs,stmp);
            ++nrows;
        }
        ifs.close();
    }

    this->resize(nrows, ncols);
    std::ifstream ifs(filename.c_str());

    if (!ifs.good()) return false;
    for (int i = 0; i < nrows; ++i) for (int j = 0; j < ncols; ++j) ifs  >> this->operator()(i,j);
    ifs.close();
    return true;
}

