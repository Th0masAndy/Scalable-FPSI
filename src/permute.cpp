#include "permute.h"

void genPermutation(u32 size, std::vector<u32> &pi)
{
    pi.resize(size);
    for (size_t i = 0; i < pi.size(); ++i) {
        pi[i] = i;
    }
    std::shuffle(pi.begin(), pi.end(), global_built_in_prg2);
    return;
}

void permute(std::vector<u32> pi, std::vector<block> &data)
{
    std::vector<block> res(data.size());
    for (size_t i = 0; i < pi.size(); ++i) {
        res[i] = data[pi[i]];
    }
    data.assign(res.begin(), res.end());
}

void senderROT(u32 numOTs, PRNG &prng, BitVector &send0, BitVector &send1, Socket &chl, u32 numThreads)
{
    // 协议中的receiver在ROT中扮演sender
    send0.resize(numOTs);
    send1.resize(numOTs);

    SilentOtExtSender sender;
    sender.configure(numOTs, numThreads, 128);

    // sender.configure(numOTs, 40, 1, SilentSecType::SemiHonest);

    std::vector<std::array<block, 2>> sendMsg(numOTs); // 用来接受sender的输出block

    // coproto::sync_wait(sender.genBaseCors(std::nullopt, prng, chl, true));//进行baseOT
    coproto::sync_wait(sender.genBaseOts(prng, chl)); // 进行baseOT

    coproto::sync_wait(sender.silentSend(sendMsg, prng, chl)); // 进行silentOT

    for (u32 i = 0; i < numOTs; i++) {
        const unsigned char *bytes = sendMsg[i][0].data();
        bool bit = (bytes[0] >> 0) & 1;
        send0[i] = bit;

        bytes = sendMsg[i][1].data();
        bit = (bytes[0] >> 0) & 1;
        send1[i] = bit;
    }
    BitVector d;
    d.resize(numOTs);

    coproto::sync_wait(chl.recv(d));
    for (u32 i = 0; i < numOTs; i++) {
        if (d[i] == 1) {
            send0[i] ^= send1[i];
            send1[i] ^= send0[i];
            send0[i] ^= send1[i];
        }
    }
}

void receiverROT(u32 numOTs, PRNG &prng, BitVector &bitV, BitVector &recv, Socket &chl, u32 numThreads)
{
    // 协议中的sender在ROT中扮演receiver
    bitV.resize(numOTs);
    recv.resize(numOTs);
    BitVector bitV0 = bitV;
    SilentOtExtReceiver receiver;
    receiver.configure(numOTs, numThreads, 128);
    // receiver.configure(numOTs, 40, 1, SilentSecType::SemiHonest);

    std::vector<block> recvMsg(numOTs); // 用来接受receiver的输出block

    // coproto::sync_wait(receiver.genBaseCors(std::nullopt, prng, chl, true));//进行baseOT
    coproto::sync_wait(receiver.genBaseOts(prng, chl));                                   // 进行baseOT
    coproto::sync_wait(receiver.silentReceive(bitV, recvMsg, prng, chl, OTType::Random)); // 进行silentOT

    for (u32 i = 0; i < numOTs; i++) {
        const unsigned char *bytes = recvMsg[i].data();
        bool bit = (bytes[0] >> 0) & 1;
        recv[i] = bit;
    }
    BitVector diff = bitV ^ bitV0;
    coproto::sync_wait(chl.send(diff));
    bitV = bitV0;
}

void perm_cons(u8 nn, std::vector<u32> &in, std::vector<u32> pi, BitVector &switch_bit)
{
    u32 elem_num = 1 << nn;
    u32 switch_num = 1 << (nn - 1);
    std::vector<u32> out;
    out.resize(elem_num);
    for (u32 i = 0; i < elem_num; i++) {
        out[i] = pi[i];
    }
    // std::vector<uint8_t> switch_bit(switch_num);//以
    switch_bit.resize((2 * nn - 1) * (1 << (nn - 1)));
    u32 switchnum_inner;
    u32 j_num;
    u32 jloop;
    std::vector<u32> all_l(elem_num);
    std::vector<u32> all_r(elem_num);

    std::vector<u32> in1(elem_num);
    std::vector<u32> out1(elem_num);

    u32 midpoint; // midpoint代表了考虑节点的位置，取值就是0，1，，，n这样的
    // nn是元素个数的对数，pi代表了这个置换，switch_bit代表了每一个位置的交换位
    for (u8 i = 0; i < nn; i++) {
        switchnum_inner = (1 << (nn - 1 - i));
        j_num = (1 << (nn - i));
        jloop = (1 << i);
        std::vector<bool> allocated(1 << nn);

        // midpoint=in[0];
        // printf("当前i=%d，此时midpoint是：%d\n",static_cast<int>(i),midpoint);
        // 需要把all_l和all_r都计算出来
        if (i == 0) { // 第一次，all用默认的in和out初始化
            for (int m = 0; m < (elem_num); m++) {
                all_r[out[m]] = m; // all_r揭示了右外侧元素的位置
                all_l[in[m]] = m;  // all_l揭示了左外侧元素的位置
            }

        } else {
            // for (int q = 0; q < (elem_num); q++)
            // {

            //     in[q]=in1[q];
            //     out[q]=out1[q];
            // }
            memcpy(in.data(), in1.data(), (elem_num) * sizeof(int));
            memcpy(out.data(), out1.data(), (elem_num) * sizeof(int));

            // 随后的每一次，都用上一次计算出来的in1和out1来更新all_l和all_r
            for (int m = 0; m < (elem_num); m++) {
                all_r[out[m]] = m;
                all_l[in[m]] = m;
            }
        }
        for (int j = 0; j < (jloop); j++) {
            if (i == (nn - 1)) { // 最后一层，需要有别的switch确定方法，直接根据是否相等即可判断switch_bit
                if (in[2 * j] == out[2 * j]) {
                    switch_bit[(nn - 1) * (switch_num) + j] = 0; // 这里的值需要改变*************************************
                } else {
                    switch_bit[(nn - 1) * (switch_num) + j] = 1;
                }

            } else {
                bool dump = 0;
                // std::vector<bool> allocated(switchnum_inner);
                // std::vector<bool> allocated(1 << nn);

                for (int k = 0; k < (switchnum_inner); k++) {
                    if (k == 0) { // 第一次嘛，直接设置第一个switch默认为0，off
                        midpoint = in[j * (j_num)];
                        allocated[midpoint] = 1;
                        switch_bit[i * (switch_num) + j * (switchnum_inner)] = 0; // 第一个switch默认为off，设置为0
                        in1[all_l[midpoint]] = midpoint;
                        in1[all_l[midpoint] + (switchnum_inner)] = in[all_l[midpoint] + 1]; // 这里刚刚改过
                        allocated[in[all_l[midpoint] + 1]] = 1;
                        // 这里
                        if (all_r[midpoint] % 2 == 1) {
                            switch_bit[(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                1; //+j*(1 << (nn-1-i))
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                            midpoint = out[all_r[midpoint] - 1];
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;

                        } else {
                            switch_bit[(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                0; // 这里
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                            midpoint = out[all_r[midpoint] + 1];
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                        }

                    } else {
                        while (allocated[midpoint] == 1) // 这段代码应该写在右侧结束后，而不是左侧开始前
                        {
                            midpoint = in[j * (j_num) + ((all_l[midpoint] - j * (j_num)) + 2) % (j_num)];
                            dump = 1;
                        }

                        if (dump == 1) { // 发生了无意义的循环，需要修正
                            dump = 0;    // 重新把循环的开关关上
                            allocated[midpoint] = 1;
                            switch_bit[i * (switch_num) + j * (switchnum_inner) + ((all_l[midpoint] % (j_num)) / 2)] = 0;
                            in1[j * (j_num) + (all_l[midpoint] - j * (j_num)) / 2 + (switchnum_inner)] = midpoint;
                            midpoint = in[all_l[midpoint] - 1]; // 这里的逻辑
                            in1[j * (j_num) + (all_l[midpoint] - j * (j_num)) / 2] = midpoint;
                            allocated[midpoint] = 1;

                            if (all_r[midpoint] % 2 == 1) {
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                        1; //+j*(1 << (nn-1-i))
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                                midpoint = out[all_r[midpoint] - 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;

                            } else {
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                        0; // 这里
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                                midpoint = out[all_r[midpoint] + 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                            }

                        } else {
                            if (all_l[midpoint] % 2 == 1) {
                                switch_bit[i * (switch_num) + j * (switchnum_inner) + ((all_l[midpoint] - j * (j_num)) / 2)] = 0;
                                in1[all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint; // zheli
                                allocated[midpoint] = 1;
                                midpoint = in[all_l[midpoint] - 1];
                                allocated[midpoint] = 1;
                                in1[all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2)] = midpoint;

                            } else {
                                switch_bit[i * (switch_num) + j * (switchnum_inner) + ((all_l[midpoint] - j * (j_num)) / 2)] = 1;
                                in1[all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                                allocated[midpoint] = 1;
                                midpoint = in[all_l[midpoint] + 1];
                                allocated[midpoint] = 1;
                                in1[(all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2))] = midpoint; // 这里有一个r
                            }

                            if (all_r[midpoint] % 2 == 0) { // zheli
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] = 0;
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2)] = midpoint;
                                midpoint = out[all_r[midpoint] + 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;

                            } else {
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] = 1;
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2)] = midpoint;
                                midpoint = out[all_r[midpoint] - 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                            }
                        }
                    }
                }
            }
        }
        // for (int q = 0; q < (elem_num); q++)
        // {
        //     // in[q]=in1[q];
        //     // out[q]=out1[q];
        // }
    }
}

void Opermute(u32 idx, BitVector &input, BitVector &output, std::vector<u32> &pi, Socket &chl, u32 numThreads)
{
    // 汇总上述函数，进行置换网络的构建，双方交互的过程
    u32 num = input.size();                   // 这是置换有用的元素的个数
    u8 nn = std::bit_width(input.size() - 1); // nn是置换网络的参数
    u32 elem_num = 1 << nn;
    u32 switch_num = (2 * nn - 1) * (1 << (nn - 1)); // 总的switch数量
    input.resize(elem_num);                          // input和output都扩展到2^nn个元素
    // std::cout<<"forrmal input:"<<idx<<" :"<<input<<std::endl;

    output.resize(elem_num);
    BitVector input0;
    input0.resize(elem_num);
    BitVector d(switch_num);
    if (idx == 0) { // sender
        // 生成置换

        std::vector<u32> norm; // 原始的顺序，0，1，2，，，
        BitVector reality;     // 记录了哪些元素是有用的bit
        norm.resize(elem_num);
        reality.resize(elem_num);
        pi.resize(elem_num);
        for (u32 i = 0; i < elem_num; i++) { // 置换的输入输出
            norm[i] = i;
            pi[i] = i;
        }

        std::shuffle(pi.begin(), pi.end(), global_built_in_prg2); // 生成随机置换
        // std::cout<<"pi0:"<<std::endl;
        // for(u32 i =0;i<elem_num;i++){
        //     std::cout<<pi[i]<<std::endl;

        // }

        for (u32 i = 0; i < elem_num; i++) { // 给reality赋值
            if (pi[i] < num) {
                reality[i] = 1;
            }
        }
        // std::cout<<"Sreality:"<<reality<<std::endl;

        // 构建网络
        BitVector network_choice;
        Timer timer1;
        timer1.setTimePoint("start");
        perm_cons(nn, norm, pi, network_choice); // 构建网络，得到每一个switch的选择位
        timer1.setTimePoint("ours");
        std::cout << timer1 << std::endl;

        // std::cout<<"pi置换网络:"<<network_choice<<std::endl;
        // std::cout<<"pi1:"<<std::endl;
        // for(u32 i =0;i<elem_num;i++){
        //     std::cout<<pi[i]<<std::endl;

        // }

        // 进行每一个switch的ROT
        bool isSender = (idx == 1); // S是OT的receiver
        PRNG prng(sysRandomSeed());
        BitVector bitV, recv;
        bitV = network_choice;
        // std::cout<<"sender的choicebit:"<<bitV<<std::endl;

        receiverROT(switch_num, prng, bitV, recv, chl, numThreads);
        // std::cout<<"sender的rc:"<<recv<<std::endl;

        // 本地计算最后的output
        coproto::sync_wait(chl.recv(d));
        // std::cout<<"sender的d:"<<d<<std::endl;

        // std::cout<<"S111"<<std::endl;

        for (u8 i = 0; i < (nn - 1); i++) { // 左到中间
            u32 group_num = 1 << i;
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (nn - 1 - i);
                u32 elem_num0 = 1 << (nn - i);
                for (u32 k = 0; k < switch_num0; k++) {
                    // std::cout<<"这是多少"<<(i*(elem_num/2)+j*switch_num0+k)<<std::endl;

                    if (bitV[i * (elem_num / 2) + j * switch_num0 + k] == 0) {
                        input0[j * elem_num0 + k] = input[j * elem_num0 + 2 * k] ^ recv[i * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + k + switch_num0] =
                            input[j * elem_num0 + 2 * k + 1] ^ recv[i * (elem_num / 2) + j * switch_num0 + k] ^ d[i * (elem_num / 2) + j * switch_num0 + k];
                    } else if (bitV[i * (elem_num / 2) + j * switch_num0 + k] == 1) {
                        input0[j * elem_num0 + k] =
                            input[j * elem_num0 + 2 * k + 1] ^ recv[i * (elem_num / 2) + j * switch_num0 + k] ^ d[i * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + k + switch_num0] = input[j * elem_num0 + 2 * k] ^ recv[i * (elem_num / 2) + j * switch_num0 + k];
                        // std::cout<<"       1209"<<std::endl;
                    }
                }
            }
            input = input0;
            // std::cout<<"Sinput状态："<<input<<std::endl;
        }
        // std::cout<<"Sinput0:"<<input<<std::endl;

        for (u32 i = 0; i < (elem_num / 2); i++) { // 中间
            // std::cout<<"i:"<<i<<std::endl;

            if (bitV[(nn - 1) * (elem_num / 2) + i] == 0) {
                input0[(i << 1)] = input[i << 1] ^ recv[(nn - 1) * (elem_num / 2) + i];
                input0[(i << 1) + 1] = input[(i << 1) + 1] ^ recv[(nn - 1) * (elem_num / 2) + i] ^ d[(nn - 1) * (elem_num / 2) + i];
            } else if (bitV[(nn - 1) * (elem_num / 2) + i] == 1) {
                input0[(i << 1)] = input[(i << 1) + 1] ^ recv[(nn - 1) * (elem_num / 2) + i] ^ d[(nn - 1) * (elem_num / 2) + i];
                input0[(i << 1) + 1] = input[i << 1] ^ recv[(nn - 1) * (elem_num / 2) + i];
            }
        }
        input = input0;
        // std::cout<<"Sinput状态："<<input<<std::endl;

        // std::cout<<"Sinput1:"<<input<<std::endl;

        for (u8 i = 0; i < (nn - 1); i++) { // 中间到右
            u32 group_num = 1 << (nn - 2 - i);
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (i + 1);
                u32 elem_num0 = 1 << (i + 2);
                for (u32 k = 0; k < switch_num0; k++) {
                    if (bitV[(nn + i) * (elem_num / 2) + j * switch_num0 + k] == 0) {
                        input0[j * elem_num0 + 2 * k] = input[j * elem_num0 + k] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + 2 * k + 1] = input[j * elem_num0 + k + switch_num0] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k] ^
                                                            d[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    } else if (bitV[(nn + i) * (elem_num / 2) + j * switch_num0 + k] == 1) {
                        input0[j * elem_num0 + 2 * k] = input[j * elem_num0 + k + switch_num0] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k] ^
                                                        d[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + 2 * k + 1] = input[j * elem_num0 + k] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    }
                }
            }
            input = input0;
            // std::cout<<"Sinput状态："<<input<<std::endl;
        }

        // output=input;
        coproto::sync_wait(chl.send(reality));

        // std::cout<<"senderoutput:"<<input<<std::endl;
        for (u32 jelly = 0, i = 0; i < elem_num; i++) {
            if (reality[i] == 1) {
                output[jelly] = input[i];
                jelly++;
            }
        }
        output.resize(num);

        std::vector<u32> pi0(num);
        for (u32 i = 0, j = 0; i < elem_num; i++) { // 对pi进行修改，由原来的numBin大小，改为numelements大小
            if (pi[i] < num) {
                pi0[j] = pi[i];
                j++;
            }
        }

        pi.resize(num);

        for (u32 i = 0; i < num; i++) {
            pi[i] = pi0[i];
        }

    } else if (idx == 1) {
        // 进行每一个switch的ROT
        bool isSender = (idx == 1); // R是OT的sender
        PRNG prng(sysRandomSeed());
        BitVector send0, send1;
        senderROT(switch_num, prng, send0, send1, chl, numThreads);

        // 本地计算最后的output
        for (u8 i = 0; i < (nn - 1); i++) { // 左到中间
            u32 group_num = 1 << i;
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (nn - 1 - i);
                u32 elem_num0 = 1 << (nn - i);
                for (u32 k = 0; k < switch_num0; k++) {
                    input0[j * elem_num0 + k] = input[j * elem_num0 + 2 * k] ^ send0[i * (elem_num / 2) + j * switch_num0 + k];
                    input0[j * elem_num0 + k + switch_num0] = input[j * elem_num0 + 2 * k] ^ send1[i * (elem_num / 2) + j * switch_num0 + k];
                    d[i * (elem_num / 2) + j * switch_num0 + k] =
                        input[j * elem_num0 + 2 * k] ^ input[j * elem_num0 + 2 * k + 1] ^ input0[j * elem_num0 + k] ^ input0[j * elem_num0 + k + switch_num0];
                } // 写到了这里
            }
            input = input0;
            // std::cout<<"Rinput状态："<<input<<std::endl;
        }

        // std::cout<<"Rinput0:"<<input<<std::endl;

        for (u32 i = 0; i < (elem_num / 2); i++) { // 中间
            input0[i << 1] = input[i << 1] ^ send0[(nn - 1) * (elem_num / 2) + i];
            input0[(i << 1) + 1] = input[i << 1] ^ send1[(nn - 1) * (elem_num / 2) + i];
            d[(nn - 1) * (elem_num / 2) + i] = input0[i << 1] ^ input0[(i << 1) + 1] ^ input[i << 1] ^ input[(i << 1) + 1];
        }
        input = input0;

        for (u8 i = 0; i < (nn - 1); i++) {    // 中间到右
            u32 group_num = 1 << (nn - 2 - i); //
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (i + 1); //
                u32 elem_num0 = 1 << (i + 2);   //
                for (u32 k = 0; k < switch_num0; k++) {
                    input0[j * elem_num0 + 2 * k] = input[j * elem_num0 + k] ^ send0[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    input0[j * elem_num0 + 2 * k + 1] = input[j * elem_num0 + k] ^ send1[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    d[(nn + i) * (elem_num / 2) + j * switch_num0 + k] =
                        input0[j * elem_num0 + 2 * k] ^ input0[j * elem_num0 + 2 * k + 1] ^ input[j * elem_num0 + k] ^ input[j * elem_num0 + k + switch_num0];
                }
            }
            input = input0;
            // std::cout<<"Rinput状态："<<input<<std::endl;
        }

        // std::cout<<"receiver的d:"<<d<<std::endl;

        coproto::sync_wait(chl.send(d));
        BitVector reality;
        reality.resize(elem_num);
        coproto::sync_wait(chl.recv(reality));
        // std::cout<<"Routput:"<<input<<std::endl;
        // std::cout<<"Rreality:"<<reality<<std::endl;

        for (u32 jelly = 0, i = 0; i < elem_num; i++) {
            if (reality[i] == 1) {
                output[jelly] = input[i];
                jelly++;
            }
        }
        output.resize(num);

        // output=input;
    }
}

void Opermute1(u32 idx, BitVector &input, BitVector &output, std::vector<u32> &pi, Socket &chl, u32 numThreads)
{
    u32 num = input.size();                     // 这是置换的元素的个数
    u8 nn = std::bit_width(input.size() - 1);   // nn是元素个数的对数上界
    u32 layers = 2 * nn - 1;                    // 一共的网络层数
    u32 switch_num_layer = num / 2;             // 每一层的switch数量
    u32 switch_num = layers * switch_num_layer; // 总的switch数量
    BitVector d(switch_num);                    // sender和recver通信传递的d值

    pi.resize(num);
    output.resize(num);
    for (u32 i = 0; i < num; i++) { // 置换的输入输出
        pi[i] = i;
    }

    std::shuffle(pi.begin(), pi.end(), global_built_in_prg2); // 生成随机置换

    std::vector<int32_t> pi1(pi.begin(), pi.end());
    Timer timer2;
    timer2.setTimePoint("start");
    WaksmanNetwork network(pi1); // 根据随机置换生成网络，sender的用来做路由，recver的用来指示，-1,0,1,2
    auto switches = network.get_waksman_network();
    timer2.setTimePoint("alibaba");
    // std::cout << timer2 << std::endl;

    if (idx == 0) {
        // 进行ROT交互

        bool isSender = (idx == 1); // S是OT的receiver
        PRNG prng(sysRandomSeed());
        BitVector recv;
        BitVector bitV(switch_num);
        // bitV.resize(switch_num);
        for (u32 i = 0; i < layers; i++) {
            for (u32 j = 0; j < switch_num_layer; j++) {
                if (switches[i][j] == 1) {
                    bitV[i * switch_num_layer + j] = 1;
                } else {
                    bitV[i * switch_num_layer + j] = 0;
                }
            }
        }

        // bitV=network_choice;
        // std::cout<<"sender的choicebit:"<<bitV<<std::endl;
        receiverROT(switch_num, prng, bitV, recv, chl, numThreads);
        // std::cout<<"sender的rc:"<<recv<<std::endl;

        coproto::sync_wait(chl.recv(d));
        // 进行本地的计算

        u32 loog = 1;
        u32 start = 0;
        std::vector<u32> group_size(switch_num_layer);
        std::vector<u32> switch_start(switch_num_layer);
        BitVector input0(num);

        // group_size.resize(1);
        // switch_start.resize(1);
        group_size[0] = num;
        switch_start[0] = 0;

        // 进行前 2logN-2 层的计算
        for (u8 i = 0; i < nn - 2; i++) { // 一层
            start = 0;
            u32 iswitchnum = i * switch_num_layer;
            for (u32 j = 0; j < loog; j++) { // 一个层里面的一个子网络
                u32 half_groupsize = group_size[j] / 2;
                for (u32 k = 0; k < half_groupsize; k++) { // 子网络中的元素计算

                    if (bitV[iswitchnum + switch_start[j] + k] == 1) {
                        input0[start + k] = input[start + k * 2 + 1] ^ recv[iswitchnum + switch_start[j] + k] ^ d[iswitchnum + switch_start[j] + k];
                        input0[start + k + half_groupsize] = input[start + k * 2] ^ recv[iswitchnum + switch_start[j] + k];
                    } else {
                        input0[start + k] = input[start + k * 2] ^ recv[iswitchnum + switch_start[j] + k];
                        input0[start + k + half_groupsize] =
                            input[start + k * 2 + 1] ^ recv[iswitchnum + switch_start[j] + k] ^ d[iswitchnum + switch_start[j] + k];
                    }
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }
                start += group_size[j];
            }
            // group_size.resize(loog<<1);
            // switch_start.resize(loog<<1);
            for (int ite = loog - 1; ite >= 0; ite--) {
                // 写入groupsize
                group_size[(ite << 1) + 1] = (group_size[ite] + 1) / 2;                        // 除2上界
                group_size[ite << 1] = group_size[ite] / 2;                                    // 除2下界
                switch_start[(ite << 1) + 1] = switch_start[ite] + (group_size[ite << 1] / 2); // 右子网络的switch起始位置
                switch_start[ite << 1] = switch_start[ite];                                    // 左子网络的switch起始位置
            }

            loog = loog << 1;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input0, input);
        }
        start = 0;

        // cout<<"wuwuwu"<<endl;
        // 进行中间三层的计算
        for (u32 i = 0; i < group_size.size(); i++) {
            u32 nn_2layer = switch_num_layer * (nn - 2);
            u32 nn_1layer = switch_num_layer * (nn - 1);
            u32 nn_0layer = switch_num_layer * nn;

            if (group_size[i] == 4) {
                // 第一层
                if (bitV[nn_2layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                    input0[start + 2] = input[start] ^ recv[nn_2layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_2layer + switch_start[i]];
                    input0[start + 2] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                }

                if (bitV[nn_2layer + switch_start[i] + 1] == 1) {
                    input0[start + 1] = input[start + 3] ^ recv[nn_2layer + switch_start[i] + 1] ^ d[nn_2layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 2] ^ recv[nn_2layer + switch_start[i] + 1];
                } else {
                    input0[start + 1] = input[start + 2] ^ recv[nn_2layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 3] ^ recv[nn_2layer + switch_start[i] + 1] ^ d[nn_2layer + switch_start[i] + 1];
                }
                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];

                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第二层
                if (bitV[nn_1layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_1layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                }

                if (bitV[nn_1layer + switch_start[i] + 1] == 1) {
                    input0[start + 2] = input[start + 3] ^ recv[nn_1layer + switch_start[i] + 1] ^ d[nn_1layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 2] ^ recv[nn_1layer + switch_start[i] + 1];
                } else {
                    input0[start + 2] = input[start + 2] ^ recv[nn_1layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 3] ^ recv[nn_1layer + switch_start[i] + 1] ^ d[nn_1layer + switch_start[i] + 1];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第三层
                if (bitV[nn_0layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 2] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_0layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start + 2] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                }

                input0[start + 2] = input[start + 1];
                input0[start + 3] = input[start + 3]; // 一般来说，这个switch可以忽略，采用直连的策略

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 3) {
                // 第一层
                if (bitV[nn_2layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_2layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_2layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第二层
                if (bitV[nn_1layer + switch_start[i]] == 1) {
                    input0[start + 1] = input[start + 2] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                    input0[start + 2] = input[start + 1] ^ recv[nn_1layer + switch_start[i]];
                } else {
                    input0[start + 1] = input[start + 1] ^ recv[nn_1layer + switch_start[i]];
                    input0[start + 2] = input[start + 2] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第三层
                if (bitV[nn_0layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_0layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 2) { // 直接对第三层的input值进行修改，但是d还是要修改到中间的层
                if (bitV[nn_1layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_1layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                // input   = std::move(input0);
                // std::swap(input, input0);
            }
            start += group_size[i];
        }
        // cout<<"wuwuwu"<<endl;

        // 进行后 2logN-2 层的计算
        loog = loog / 2;
        for (u8 i = 0; i < nn - 2; i++) { // 后半部分网络的计算
            for (int ite = 0; ite < loog; ite++) {
                // 写入groupsize
                group_size[ite] = group_size[ite << 1] + group_size[(ite << 1) + 1]; // 除2上界
                switch_start[ite] = switch_start[ite << 1];                          // 左子网络的switch起始位置
            }
            // group_size.resize(loog);
            // switch_start.resize(loog);

            start = 0;
            u32 nn_1ilayer = (nn + 1 + i) * switch_num_layer;
            for (u32 j = 0; j < loog; j++) {
                for (u32 k = 0; k < group_size[j] / 2; k++) { // 子网络中的元素计算
                    if (bitV[nn_1ilayer + switch_start[j] + k] == 1) {
                        input0[start + 2 * k] =
                            input[start + k + (group_size[j] / 2)] ^ recv[nn_1ilayer + switch_start[j] + k] ^ d[nn_1ilayer + switch_start[j] + k];
                        input0[start + 2 * k + 1] = input[start + k] ^ recv[nn_1ilayer + switch_start[j] + k];
                    } else {
                        input0[start + 2 * k] = input[start + k] ^ recv[nn_1ilayer + switch_start[j] + k];
                        input0[start + 2 * k + 1] =
                            input[start + k + (group_size[j] / 2)] ^ recv[nn_1ilayer + switch_start[j] + k] ^ d[nn_1ilayer + switch_start[j] + k];
                    }
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }
                start += group_size[j];
            }
            loog = loog / 2;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input, input0);
        }
        output = input;

    } else if (idx == 1) {
        // 进行ROT交互
        bool isSender = (idx == 1); // R是OT的sender
        PRNG prng(sysRandomSeed());
        BitVector send0, send1;
        senderROT(switch_num, prng, send0, send1, chl, numThreads);

        // 进行本地的计算

        // 进行前 2logN-2 层的计算
        u32 loog = 1;
        u32 start = 0;
        std::vector<u32> group_size(switch_num_layer);
        std::vector<u32> switch_start(switch_num_layer);
        BitVector input0(num);

        // group_size.resize(1);
        // switch_start.resize(1);
        group_size[0] = num;
        switch_start[0] = 0;
        // cout<<"wuwuwu"<<endl;

        for (u8 i = 0; i < nn - 2; i++) { // 一层
            start = 0;

            for (u32 j = 0; j < loog; j++) { // 一个层里面的一个子网络
                u32 half_group = group_size[j] / 2;
                for (u32 k = 0; k < half_group; k++) {                                                                         // 子网络中的元素计算
                    input0[start + k] = input[start + k * 2] ^ send0[switch_num_layer * i + switch_start[j] + k];              // u0 = x0 − r0,
                    input0[start + k + half_group] = input[start + k * 2] ^ send1[switch_num_layer * i + switch_start[j] + k]; // u1 = x0 − r1
                    d[switch_num_layer * i + switch_start[j] + k] = input[start + k * 2] ^ input[start + k * 2 + 1] ^
                                                                    send0[switch_num_layer * i + switch_start[j] + k] ^
                                                                    send1[switch_num_layer * i + switch_start[j] + k]; // d = (x0 + x1) − (u0 + u1)
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }
                start += group_size[j];
            }
            // group_size.resize(loog<<1);
            // switch_start.resize(loog<<1);
            for (int ite = loog - 1; ite >= 0; ite--) {
                // 写入groupsize
                group_size[(ite << 1) + 1] = (group_size[ite] + 1) / 2;                        // 除2上界
                group_size[ite << 1] = group_size[ite] / 2;                                    // 除2下界
                switch_start[(ite << 1) + 1] = switch_start[ite] + (group_size[ite << 1] / 2); // 右子网络的switch起始位置
                switch_start[ite << 1] = switch_start[ite];                                    // 左子网络的switch起始位置
            }
            loog = loog << 1;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input, input0);
        }
        start = 0;
        // cout<<"wuwuwu"<<endl;

        // 进行中间三层的计算,直接去掉外轮递归，尽可能更详细的写这部分
        for (u32 i = 0; i < group_size.size(); i++) {
            u32 nn_2layer = switch_num_layer * (nn - 2);
            u32 nn_1layer = switch_num_layer * (nn - 1);
            u32 nn_0layer = switch_num_layer * nn;
            if (group_size[i] == 4) {
                // 第一层
                input0[start] = input[start] ^ send0[nn_2layer + switch_start[i]];
                input0[start + 1] = input[start + 2] ^ send0[nn_2layer + switch_start[i] + 1];
                input0[start + 2] = input[start] ^ send1[nn_2layer + switch_start[i]];
                input0[start + 3] = input[start + 2] ^ send1[nn_2layer + switch_start[i] + 1];
                d[nn_2layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_2layer + switch_start[i]] ^ send1[nn_2layer + switch_start[i]];
                d[nn_2layer + switch_start[i] + 1] =
                    input[start + 2] ^ input[start + 3] ^ send0[nn_2layer + switch_start[i] + 1] ^ send1[nn_2layer + switch_start[i] + 1];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第二层
                input0[start] = input[start] ^ send0[nn_1layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_1layer + switch_start[i]];
                d[nn_1layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_1layer + switch_start[i]] ^ send1[nn_1layer + switch_start[i]];

                input0[start + 2] = input[start + 2] ^ send0[nn_1layer + switch_start[i] + 1];
                input0[start + 3] = input[start + 2] ^ send1[nn_1layer + switch_start[i] + 1];
                d[nn_1layer + switch_start[i] + 1] =
                    input[start + 2] ^ input[start + 3] ^ send0[nn_1layer + switch_start[i] + 1] ^ send1[nn_1layer + switch_start[i] + 1];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第三层
                input0[start] = input[start] ^ send0[nn_0layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_0layer + switch_start[i]];
                d[nn_0layer + switch_start[i]] = input[start] ^ input[start + 2] ^ send0[nn_0layer + switch_start[i]] ^ send1[nn_0layer + switch_start[i]];

                input0[start + 2] = input[start + 1];
                input0[start + 3] = input[start + 3]; // 一般来说，这个switch可以忽略，采用直连的策略
                // d[switch_num_layer*nn+switch_start[i]+1]=input[start+2]^input[start+3]^send0[switch_num_layer*nn+switch_start[i]+1]^send1[switch_num_layer*nn+switch_start[i]+1];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 3) {
                // 第一层
                input0[start] = input[start] ^ send0[nn_2layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_2layer + switch_start[i]];
                d[nn_2layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_2layer + switch_start[i]] ^ send1[nn_2layer + switch_start[i]];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第二层
                input0[start + 1] = input[start + 1] ^ send0[nn_1layer + switch_start[i]];
                input0[start + 2] = input[start + 1] ^ send1[nn_1layer + switch_start[i]];
                d[nn_1layer + switch_start[i]] = input[start + 1] ^ input[start + 2] ^ send0[nn_1layer + switch_start[i]] ^ send1[nn_1layer + switch_start[i]];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // 第三层
                input0[start] = input[start] ^ send0[nn_0layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_0layer + switch_start[i]];
                d[nn_0layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_0layer + switch_start[i]] ^ send1[nn_0layer + switch_start[i]];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 2) { // 直接对第三层的input值进行修改，但是d还是要修改到中间的层
                input0[start] = input[start] ^ send0[nn_1layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_1layer + switch_start[i]];
                d[nn_1layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_1layer + switch_start[i]] ^ send1[nn_1layer + switch_start[i]];

                input = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                // input   = std::move(input0);
                // std::swap(input, input0);
            }
            start += group_size[i];
        }
        // cout<<"wuwuwu"<<endl;

        // 进行后 2logN-2 层的计算
        loog = loog / 2;

        for (u8 i = 0; i < nn - 2; i++) {
            u32 nn_1iswitchnum = (nn + 1 + i) * switch_num_layer;

            for (int ite = 0; ite < loog; ite++) {
                // 写入groupsize
                group_size[ite] = group_size[ite << 1] + group_size[(ite << 1) + 1];
                switch_start[ite] = switch_start[ite << 1];
            }
            // group_size.resize(loog);
            // switch_start.resize(loog);

            start = 0;
            for (u32 j = 0; j < loog; j++) { // 一个层里面的一个子网络

                for (u32 k = 0; k < group_size[j] / 2; k++) {                                               // 子网络中的元素计算
                    input0[start + 2 * k] = input[start + k] ^ send0[nn_1iswitchnum + switch_start[j] + k]; // u0 = x0 − r0,
                    input0[start + 2 * k + 1] = input[start + k] ^ send1[nn_1iswitchnum + switch_start[j] + k];
                    d[nn_1iswitchnum + switch_start[j] + k] = input[start + k] ^ input[start + k + (group_size[j] / 2)] ^
                                                              send0[nn_1iswitchnum + switch_start[j] + k] ^ send1[nn_1iswitchnum + switch_start[j] + k];
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }

                start += group_size[j];
            }
            loog = loog / 2;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input, input0);
        }
        output = input;
        // 最后的通信阶段
        coproto::sync_wait(chl.send(d));
    }
}