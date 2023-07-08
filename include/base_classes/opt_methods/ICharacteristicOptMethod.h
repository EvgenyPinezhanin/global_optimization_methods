#ifndef I_CHARACTERISTIC_OPT_METHOD_H_
#define I_CHARACTERISTIC_OPT_METHOD_H_

class ICharacteristicOptMethod {
protected:
    double t;

    virtual void calcCharacteristic() = 0;

public:
    ICharacteristicOptMethod() : t(0) {};
};

#endif // I_CHARACTERISTIC_OPT_METHOD_H_