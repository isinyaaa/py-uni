double[2] skck(double a, double b)
{
    double dem = pow(a*a + b*b, 1/2);

    double skck[2];
    skck[0] = a/dem;
    skck[1] = -b/dem;

    return skck;
}


double [2] givens_rotation
(double alphas[], double betas[], double gammas[], double deltas[], int betas_size, int i)
{
    int j = i + 1;

    double skck = skck(alphas[i], betas[i]);
    double sk = skck[0];
    double ck = skck[1];

    double prev_alpha, prev_beta, prev_gamma, prev_delta;

    for (int k = 0; k < 3; k++) {
        if (k == 2 && i != betas_size - 1) {
            prev_delta = deltas[i];
            prev_gamma = gammas[j];
            gammas[j] = sk * prev_delta + ck * prev_gamma;
            gammas[i] = ck * prev_delta - sk * prev_gamma;
        } else if (k == 1) {
            prev_gamma = gammas[i];
            prev_alpha = alphas[j];
            gammas[i] = ck * prev_gamma - sk * prev_alpha;
            alphas[j] = sk * prev_gamma + ck * prev_alpha;
        } else if (k == 0) {
            prev_alpha = alphas[i];
            prev_beta = betas[i];
            betas[i]  = sk * prev_alpha + ck * prev_beta;
            alphas[i] = ck * prev_alpha - sk * prev_beta;
        }
    }

    return {sk, ck};
}

int sign(double n)
{
    return n < 0 ? -1 : 1;
}

void RQ
(double alphas[], double betas[], double gammas[], int alpha_size, int gamma_size, double vec_ck[], double vec_sk[])
{
    int j;
    double sk, ck;
    double prev_alpha, prev_beta;

    for (int i = 0; i < alpha_size - 1; i++)
    {
        j = i + 1;
        sk = vec_sk[i];
        ck = vec_ck[i];

        prev_alpha = alphas[i];
        alphas[i] = sk * gammas[i] + ck * prev_alpha;

        prev_alpha = alphas[j];
        prev_beta = betas[i];
        betas[i] = ck * prev_beta - sk * prev_alpha;
        alphas[j] = ck * prev_alpha + sk * prev_beta;
    }

    for (int i = 0; i < gamma_size; i++)
        gammas[i] = betas[i];
}

void VQ
(double V[][], int v_size, double vec_ck[], double vec_sk)
{
    int j;
    double sk, ck;
    double prev_vki, prev_vkj;

    for (int i = 0; i < v_size - 1; i++)
    {
        j = i + 1;
        sk = vec_sk[i];
        ck = vec_ck[i];

        for (int k = 0; k < v_size; k++)
        {
            prev_vki = V[k][i];
            prev_vkj = V[k][j];
            V[k][i] = ck * prev_vkj + sk * prev_vki;
            V[k][j] = ck * prev_vki - sk * prev_vkj;
        }
    }
}

