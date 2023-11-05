#include "rid-solomon.h"

#include <string.h>

unsigned int primes[] = {
    PRIM(0), PRIM(1), PRIM(2), PRIM(3),
    PRIM(4), PRIM(5), PRIM(6), PRIM(7),
    PRIM(8), PRIM(9), PRIM(10), PRIM(11),
    PRIM(12), PRIM(13), PRIM(14), PRIM(15),
    PRIM(16), PRIM(17), PRIM(18), PRIM(19),
    PRIM(20), PRIM(21), PRIM(22), PRIM(23),
    PRIM(24), PRIM(25), PRIM(26), PRIM(27),
    PRIM(28), PRIM(29), PRIM(30), PRIM(31)
};

void FindPrimePolys(std::ostream* out, int fieldPower, int limit)
{
    GaloisField gf(fieldPower);
    int primesFound = 0;
    RS_WORD fieldCharacteristic = ((RS_WORD)1 << fieldPower) - 1, fieldCharacteristicNext = ((RS_WORD)1 << (fieldPower + 1)) - 1;
    for (RS_WORD i = fieldCharacteristic + 2; (fieldCharacteristicNext == 0 ? i > 0 : i < fieldCharacteristicNext); i += 2)
    {
        unsigned int x = 2; // пропуск первой итерации
        bool conflict = false;
        for (unsigned int j = 1; j < fieldCharacteristic; j++)
        {
            x <<= 1;
            if (x > fieldCharacteristic)
            {
                x ^= i;
            }
            if (x == 2) // код циклический, 2 всегда будет повторяться первым
            {
                conflict = true;
                break;
            }
        }
        if (!conflict)
        {
            *out << std::hex << i << std::endl;
            primesFound++;
            if (primesFound >= limit)
            {
                return;
            }
        }
    }
}

GaloisField::GaloisField(int fieldPower)
{
    _characteristic = ((RS_WORD)1 << fieldPower) - 1;
    _fieldPower = fieldPower;
    _primitivePoly = primes[fieldPower];
    // инициализация таблиц
    unsigned int val = 1;
    _powTable = (RS_WORD*)malloc(sizeof(RS_WORD) * _characteristic * 2);
    _logTable = (RS_WORD*)malloc(sizeof(RS_WORD) * (_characteristic + 1));
    _powTable[0] = val;
    _logTable[0] = 0;
    _logTable[1] = 0;
    for (RS_WORD i = 1; i < _characteristic; i++)
    {
        val <<= 1;
        if (val > _characteristic)
        {
            val ^= _primitivePoly;
        }
        _powTable[i] = (RS_WORD)val;
        _logTable[(RS_WORD)val] = i;
    }
    for (RS_WORD i = _characteristic; i < _characteristic * 2; i++)
    {
        _powTable[i] = _powTable[i - _characteristic];
    }
    /* //Debug вывод
    for (unsigned int i = 0; i < this->_characteristic * 2; i++)
        {
            std::cout << "2^" << i << "=" << _powTable[i] << std::endl;
        }
        for (unsigned int i = 0; i < this->_characteristic + 1; i++)
        {
            std::cout << "log" << i << "=" << _logTable[i] << std::endl;
        }
    */
}

GaloisField::~GaloisField()
{
    free(_powTable);
    free(_logTable);
}

RS_WORD GaloisField::multNoLUT(RS_WORD a, RS_WORD b)
{
    RS_WORD ret = 0;
    while (b > 0)
    {
        if (b & 1) //if odd
        {
            ret ^= a;
        }
        b >>= 1;
        a <<= 1;
        if (a > _characteristic)
        {
            a ^= _primitivePoly;
        }
    }
    return ret;
}

inline RS_WORD GaloisField::mult(RS_WORD a, RS_WORD b)
{
    return (a == 0 || b == 0) ? 0 : _powTable[_logTable[a] + _logTable[b]];
}
inline RS_WORD GaloisField::div(RS_WORD a, RS_WORD b)
{
    return a == 0 ? 0 : (b == 0) ? -1 : _powTable[_logTable[a] - _logTable[b] + _characteristic];
}

inline RS_WORD GaloisField::pow(RS_WORD x, RS_WORD power)
{
    return _powTable[(_logTable[x] * power) % _characteristic];
}

inline RS_WORD GaloisField::inv(RS_WORD x)
{
    return _powTable[_characteristic - _logTable[x]];
}

inline RS_WORD GaloisField::sqrt(RS_WORD x)
{
    return _logTable[x] % 2 ? _powTable[(_logTable[x] + _characteristic) / 2] : _powTable[_logTable[x] / 2];
}

Poly::Poly()
{
    this->init();
}

Poly::Poly(int n, RS_WORD* data)
{
    this->init();
    this->setCopy(n, data);
}

Poly::~Poly()
{
    if (_coef)
    {
        free(_coef);
    }
}

void Poly::init()
{
    _n = 0;
    _coef = nullptr;
}

void Poly::setCopy(int n, RS_WORD* coef)
{
    if (n > _n)
    {
        if (_coef)
        {
            free(_coef);
        }
        _coef = (RS_WORD*)malloc(sizeof(RS_WORD) * n);
    }
    _n = n;
    if (coef)
    {
        memcpy(_coef, coef, sizeof(RS_WORD) * n);
    } else
    {
        memset(_coef, 0, sizeof(RS_WORD) * n);
    }
}

void Poly::setRef(int n, RS_WORD* coef)
{
    if (_coef)
    {
        free(_coef);
    }
    _n = n;
    _coef = coef;
}

void Poly::print()
{
    std::cout << "Poly(n=" << _n << ")";
    if (_n > 0)
    {
        std::cout << "[" << std::setw(3) << std::hex << (RS_WORD)_coef[0];
    }
    for (int i = 1; i < _n; i++)
    {
        std::cout << ", " << std::setw(3) << std::hex << (RS_WORD)_coef[i];
    }
    if (_n > 0)
    {
        std::cout << "]" << std::dec << std::endl;
    }
}

Poly* Poly_Create(int n, RS_WORD* coef)
{
    Poly* poly = (Poly*)malloc(sizeof(Poly));
    poly->init();
    poly->setCopy(n, coef);
    return poly;
}

void Poly_Free(Poly* poly)
{
    free(poly->_coef);
    free(poly);
}

void Poly_Add(Poly* out, Poly* a, Poly* b)
{
    int n = MAX(a->_n, b->_n);
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * n);
    memset(temp, 0, sizeof(RS_WORD) * n);
    for (int i = 0; i < a->_n; i++)
    {
        temp[i + n - a->_n] = a->_coef[i];
    }
    for (int i = 0; i < b->_n; i++)
    {
        temp[i + n - b->_n] ^= b->_coef[i];
    }
    out->setRef(n, temp);
}

void Poly_Scale(Poly* out, Poly* in, RS_WORD scale, GaloisField* gf)
{
    if (out == in)
    {
        for (int i = 0; i < in->_n; i++)
        {
            in->_coef[i] = gf->mult(in->_coef[i], scale);
        }
    } else
    {
        RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * in->_n);
        for (int i = 0; i < in->_n; i++)
        {
            temp[i] = gf->mult(in->_coef[i], scale);
        }
        out->setRef(in->_n, temp);
    }
}

void Poly_Mult(Poly* out, Poly* a, Poly* b, GaloisField* gf)
{
    int n = a->_n + b->_n - 1;
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * n);
    memset(temp, 0, sizeof(RS_WORD) * n);
    for (int i = 0; i < a->_n; i++)
    {
        for (int j = 0; j < b->_n; j++)
        {
            temp[i + j] ^= gf->mult(a->_coef[i], b->_coef[j]);
        }
    }
    out->setRef(n, temp);
}

void Poly_Div(Poly* result, Poly* quotient, Poly* remainder, Poly* a, Poly* b, GaloisField* gf)
{
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* a->_n);
    RS_WORD normalizer = b->_coef[0];
    memcpy(temp, a->_coef, sizeof(RS_WORD) * a->_n);
    for (int i = 0; i < a->_n - b->_n + 1; i++)
    {
        temp[i] = gf->div(temp[i], normalizer);
        RS_WORD coef = temp[i];
        if (coef != 0)
        {
            for (int j = 1; j < b->_n; j++)
            {
                if (b->_coef[j] != 0)
                {
                    temp[i + j] ^= gf->mult(b->_coef[j], coef);
                }
            }
        }
    }
    if (result)
    {
        result->setCopy(a->_n, temp);
    }
    int separator = a->_n - b->_n + 1;
    if (quotient)
    {
        quotient->setCopy(separator, temp);
    }
    if (remainder)
    {
        remainder->setCopy(b->_n - 1, temp + separator);
    }
    free(temp);
}

RS_WORD Poly_Eval(Poly* poly, RS_WORD x, GaloisField* gf)
{
    RS_WORD y = poly->_coef[0];
    for (int i = 1; i < poly->_n; i++)
    {
        y = gf->mult(y, x) ^ poly->_coef[i];
    }
    return y;
}

void Poly_ChienSearch(std::vector<unsigned int>* out, Poly* poly, int max, GaloisField* gf)
{
    //this seems unnecessary because all multiplications are performed via lookup table anyway
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* poly->_n);
    memcpy(temp, poly->_coef, sizeof(RS_WORD) * poly->_n);
    for (int i = 0; i < max; i++)
    {
        RS_WORD sum = 0;
        for (int j = 0; j < poly->_n; j++)
        {
            sum ^= temp[j];
            temp[j] = gf->mult(temp[j], gf->_powTable[poly->_n - j - 1]);
        }
        if (!sum)
        {
            out->push_back(i);
        }
    }
    free(temp);
}

void Poly_Pad(Poly* poly, int left, int right)
{
    int n = poly->_n + left + right;
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* n);
    memset(temp, 0, sizeof(RS_WORD) * left);
    memcpy(temp + left, poly->_coef, sizeof(RS_WORD) * poly->_n);
    memset(temp + (left + poly->_n), 0, sizeof(RS_WORD) * right);
    poly->setRef(n, temp);
}

void Poly_Trim(Poly* poly, int left, int right)
{
    int n = poly->_n - left - right;
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* n);
    memcpy(temp, poly->_coef + left, sizeof(RS_WORD) * n);
    poly->setRef(n, temp);
}

void Poly_Append(Poly* out, Poly* a, Poly* b)
{
    int n = a->_n + b->_n;
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* n);
    memcpy(temp, a->_coef, sizeof(RS_WORD)* a->_n);
    memcpy(temp + a->_n, b->_coef, sizeof(RS_WORD)* b->_n);
    out->setRef(n, temp);
}

void Poly_Reverse(Poly* out, Poly* in)
{
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* in->_n);
    for (int i = 0; i < in->_n; i++)
    {
        temp[i] = in->_coef[in->_n - i - 1];
    }
    out->setRef(in->_n, temp);
}

ReedSolomon::ReedSolomon(int fieldPower) : _gf(fieldPower)
{

}

void ReedSolomon::createGenerator(Poly* out, int nsym)
{
    out->setCopy(1, nullptr);
    out->_coef[0] = 1;
    Poly factor(2, nullptr);
    factor._coef[0] = 1;
    for (int i = 0; i < nsym; i++)
    {
        factor._coef[1] = _gf._powTable[i];
        Poly_Mult(out, out, &factor, &_gf);
    }
}

void ReedSolomon::encode(RS_WORD* out, RS_WORD* data, int k, int nsym)
{
    Poly msg(k, data);
    Poly generator, remainder;
    this->createGenerator(&generator, nsym);
    Poly_Pad(&msg, 0, nsym);
    Poly_Div(nullptr, nullptr, &remainder, &msg, &generator, &_gf);
    memcpy(out, data, sizeof(RS_WORD) * k);
    memcpy(out + k, remainder._coef, sizeof(RS_WORD) * remainder._n);
}

void ReedSolomon::calcSyndromes(Poly* out, Poly* msg, int nsym)
{
    RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * (nsym + 1));
    for (int i = 0; i < nsym; i++)
    {
        temp[nsym - i - 1] = Poly_Eval(msg, _gf._powTable[i], &_gf);
    }
    temp[nsym] = 0; //pad
    out->setRef(nsym + 1, temp);
}

bool ReedSolomon::checkSyndromes(Poly* synd)
{
    for (int i = 0; i < synd->_n; i++)
    {
        if (synd->_coef[i])
        {
            return false;
        }
    }
    return true;
}

void ReedSolomon::findErrataLocator(Poly* out, std::vector<unsigned int>* errPos)
{
    out->setCopy(1, nullptr);
    out->_coef[0] = 1;
    Poly factor(2, nullptr);
    factor._coef[1] = 1;
    for (unsigned int i : *errPos)
    {
        factor._coef[0] = _gf._powTable[i];
        Poly_Mult(out, out, &factor, &_gf);
    }
}

void ReedSolomon::findErrorEvaluator(Poly* out, Poly* synd, Poly* errLoc, int nsym)
{
    Poly_Mult(out, synd, errLoc, &_gf); //synd lul
    Poly_Trim(out, out->_n - nsym, 0);
}

bool ReedSolomon::correctErrata(Poly* msg, Poly* synd, std::vector<unsigned int>* errPos)
{
    std::vector<unsigned int> coefPos(0);
    for (unsigned int i : *errPos)
    {
        coefPos.push_back(msg->_n - 1 - i);
    }
    Poly errLoc, errEval;
    this->findErrataLocator(&errLoc, &coefPos);
    this->findErrorEvaluator(&errEval, synd, &errLoc, errLoc._n); //TODO determine if correct
    //Poly_Reverse(errEval, errEval); //reverse it for later use
    std::vector<RS_WORD> x(coefPos.size());
    for (int i = 0; (unsigned long long)i < x.size(); i++)
    {
        x[i] = _gf._powTable[coefPos[i]]; //TODO determine if correct
    }
    Poly e(msg->_n, nullptr);
    for (int i = 0; (unsigned long long)i < x.size(); i++)
    {
        RS_WORD xi = _gf._powTable[_gf._characteristic - coefPos[i]]; //TODO determine if equivalent to GaloisField::Inv(x[i])
        RS_WORD errLocPrime = 1;
        for (int j = 0; (unsigned long long)j < x.size(); j++)
        {
            if (j != i)
            {
                errLocPrime = _gf.mult(errLocPrime, 1 ^ _gf.mult(xi, x[j]));
            }
        }
        if (errLocPrime == 0)
        {
            return false;
        }
        RS_WORD y = _gf.mult(x[i], Poly_Eval(&errEval, xi, &_gf)); //errEval is still reversed from earlier
        //TODO determine if equivalent to GaloisField::Mult(GaloisField::Pow(xi, 1), y)
        e._coef[errPos->at(i)] = _gf.div(y, errLocPrime); //magnitude
    }
    Poly_Add(msg, msg, &e);
    return true;
}

bool ReedSolomon::findErrorLocator(Poly* out, Poly* synd, int nsym, Poly* eraseLoc, int eraseCount)
{
    //this spits out a polynomial in reverse order but i dont know why
    RS_WORD init = 1;
    Poly errLoc(1, &init);
    Poly oldLoc(1, &init);
    Poly temp;
    if (eraseLoc)
    {
        errLoc.setCopy(eraseLoc->_n, eraseLoc->_coef);
        oldLoc.setCopy(eraseLoc->_n, eraseLoc->_coef);
    }
    int syndShift = 0; //TODO optimize
    for (int i = nsym - eraseCount - 1; i >= 0; i--)
    {
        int K = i + syndShift + eraseCount;
        RS_WORD delta = synd->_coef[K];
        for (int j = 1; j < errLoc._n; j++)
        {
            delta ^= _gf.mult(errLoc._coef[errLoc._n - j - 1], synd->_coef[K + j]);
        }
        Poly_Pad(&oldLoc, 0, 1); //TODO optimize
        if (delta != 0)
        {
            if (oldLoc._n > errLoc._n)
            { //TODO isn't this always the case?
                Poly_Scale(&temp, &oldLoc, delta, &_gf);
                Poly_Scale(&oldLoc, &errLoc, _gf.inv(delta), &_gf);
                errLoc.setCopy(temp._n, temp._coef);
            }
            Poly_Scale(&temp, &oldLoc, delta, &_gf);
            Poly_Add(&errLoc, &errLoc, &temp);
        }
    }
    int leading = 0;
    for (; errLoc._coef[leading] == 0; leading++);
    Poly_Trim(&errLoc, leading, 0);
    int errs = errLoc._n - 1;
    out->setCopy(errLoc._n, errLoc._coef);
    if (errs * 2 - eraseCount > nsym)
    {
        return false;
    }
    return true;
}

bool ReedSolomon::findErrors(std::vector<unsigned int>* out, Poly* errLoc, int n)
{
    int errs = errLoc->_n - 1;
    Poly revErrLoc;
    Poly_Reverse(&revErrLoc, errLoc);
    if (errLoc->_n == 1)
    {
        //do something special here? idk
    }
    else if (errLoc->_n == 2)
    { //linear equation
        out->push_back(_gf._logTable[_gf.div(errLoc->_coef[0], errLoc->_coef[1])]);
    }
    else
    {
        Poly_ChienSearch(out, &revErrLoc, n, &_gf);
    }
    if (out->size() != (unsigned long long)errs)
    {
        // Too many (or few) errors found by Chien Search for the errata locator polynomial!
        return false;
    }
    //map to string pos
    for (RS_WORD i = 0; i < out->size(); i++)
    {
        if (out->at(i) >= (unsigned long long)n) //clearly something messed up
        {
            return false;
        }
        (*out)[i] = n - out->at(i) - 1;
    }
    return true;
}

void ReedSolomon::forneySyndromes(Poly* out, Poly* synd, std::vector<unsigned int>* pos, int n)
{
    Poly fsynd(synd->_n - 1, synd->_coef);
    if (pos)
    {
        for (unsigned int i : *pos)
        {
            RS_WORD rev = (RS_WORD)n - i - 1;
            RS_WORD x = _gf._powTable[rev];
            for (int j = fsynd._n - 2; j >= 0; j--)
            {
                fsynd._coef[j + 1] = _gf.mult(fsynd._coef[j + 1], x) ^ fsynd._coef[j];
            }
        }
    }
    //fsynd.coef[fsynd.n - 1] = 0;
    out->setCopy(fsynd._n, fsynd._coef);
}

bool ReedSolomon::decode(RS_WORD* wholeOut, RS_WORD* out, RS_WORD* data, int k, int nsym, std::vector<unsigned int>* erasePos, bool debug)
{
    Poly synd;
    Poly msg(k + nsym, data);
    if (erasePos)
    {
        if (erasePos->size() > (unsigned long long)nsym)
        {
            if (debug) std::cout << "too many erasures to be corrected" << std::endl;
            return false;
        } else
        {
            for (unsigned int i : *erasePos)
            {
                msg._coef[i] = 0;
            }
        }
    }
    this->calcSyndromes(&synd, &msg, nsym);
    if (this->checkSyndromes(&synd))
    {
        if (debug) std::cout << "errors detected, locating" << std::endl;
        Poly fsynd, errLoc;
        this->forneySyndromes(&fsynd, &synd, erasePos, k + nsym);
        bool canLocate = this->findErrorLocator(&errLoc, &fsynd, nsym, nullptr, erasePos ? erasePos->size() : 0);
        if (!canLocate)
        {
            if (debug) std::cout << "too many errors to locate!" << std::endl;
            return false;
        }
        std::vector<unsigned int> pos;
        canLocate = this->findErrors(&pos, &errLoc, k + nsym);
        if (!canLocate || !(pos.size() || (erasePos && erasePos->size())))
        {
            if(debug) std::cout << "errors unable to be located!" << std::endl;
            return false;
        }
        if (debug)
        {
            if (pos.size())
            {
                std::cout << "additional errors detected at ";
                for_each(pos.begin(), pos.end(), [](unsigned int e) {std::cout << (int)e << ", "; });
            }
            std::cout << "correcting" << std::endl;
        }
        if (erasePos)
        {
            pos.insert(pos.begin(), erasePos->begin(), erasePos->end());
        }
        bool success = this->correctErrata(&msg, &synd, &pos);
        if (!success)
        {
            if(debug) std::cout << "decode failure!" << std::endl;
            return false;
        }
        if(debug) std::cout << "errors corrected" << std::endl;
    }
    if (wholeOut)
    {
        memcpy(wholeOut, msg._coef, sizeof(RS_WORD) * (k + nsym));
    }
    if (out)
    {
        memcpy(out, msg._coef, sizeof(RS_WORD) * k);
    }
    return true;
}

int code_file(FILE *input_file, FILE *coded_file) {
    RS_WORD read_buf[BUFSIZE];  //буфер чтения
    RS_WORD write_buf[BUFSIZE + NCHECKED]; //буфер записи
    memset(read_buf, '\0', BUFSIZE);
    memset(write_buf, '\0', BUFSIZE + NCHECKED);
    ReedSolomon rs(BUFSIZE);

    fseek(input_file,0,SEEK_END); //перемещение в конец файла чтобы узнать его размер
    int file_size=ftell(input_file); //узнали размер файла
    rewind(input_file); //вернули "курсор" в начало файла

    //в цикле кодируем вообщение по частям из буфера
    for (unsigned int i=0; i < file_size/BUFSIZE + (file_size%BUFSIZE ? 1 : 0) ; i++){
        fread (read_buf, BUFSIZE, 1, input_file);

        Poly r(BUFSIZE + NCHECKED, read_buf);
        rs.encode(r._coef, read_buf, BUFSIZE, NCHECKED);

        fwrite(r._coef, BUFSIZE + NCHECKED, 1, coded_file);
    }

    fclose(input_file);
    fclose(coded_file);
    std::cout << "CODED SUCCESSFULLY" << std::endl;
    return 0;
}

int decode_file(FILE *input_file, FILE *decoded_file)
{
    RS_WORD read_buf[BUFSIZE];  //буфер чтения
    RS_WORD write_buf[BUFSIZE]; //буфер записи
    memset(read_buf, '\0', BUFSIZE);
    memset(write_buf, '\0', BUFSIZE);
    ReedSolomon rs(BUFSIZE);

    fseek(input_file,0,SEEK_END); //перемещение в конец файла чтобы узнать его размер
    int file_size=ftell(input_file); //узнали размер файла
    rewind(input_file); //вернули "курсор" в начало файла

    //в цикле декодируем вообщение по частям из буфера
    for (unsigned int i=0; i < file_size/(BUFSIZE+NCHECKED) + (file_size%(BUFSIZE+NCHECKED) ? 1 : 0); i++){
        fread(read_buf, BUFSIZE + NCHECKED, 1, input_file);

        Poly e(BUFSIZE + NCHECKED, read_buf);
        Poly d(BUFSIZE, write_buf);
        bool f = rs.decode(d._coef, nullptr, e._coef, BUFSIZE, NCHECKED, nullptr, true);
        if (f) {
            fwrite(d._coef, BUFSIZE, 1, decoded_file);
        } else {
            std::cout << "ERROR DECODING" << std::endl;
        }
    }

    fclose(input_file);
    fclose(decoded_file);
    std::cout << "DECODED SUCCESSFULLY" << std::endl;
    return 0;
}

