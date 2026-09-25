import os
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Fractional Brownian surface
# ============================================================

def generate_fbm_surface(N, H, seed=None):
    """
    Generate a fractional Brownian surface using spectral synthesis.

    Parameters
    ----------
    N : int
        Linear size of the surface.
    H : float
        Hurst exponent.
    seed : int or None
        Random seed.

    Returns
    -------
    surface : ndarray
        Normalized fractional Brownian surface in the range [0, 1].
    """

    rng = np.random.default_rng(seed)

    # Frequency grid
    fx = np.fft.fftfreq(N)[:, None]
    fy = np.fft.fftfreq(N)[None, :]
    f = np.sqrt(fx**2 + fy**2)

    # Avoid division by zero at the origin
    f[0, 0] = 1.0

    # Power spectrum
    spectrum = 1.0 / f**(2 * H + 1)

    # Complex Gaussian random field
    noise = (
        rng.normal(size=(N, N))
        + 1j * rng.normal(size=(N, N))
    )

    # Spectral filtering
    F = noise * np.sqrt(spectrum)

    # Transform to spatial domain
    surface = np.fft.ifft2(F).real

    # Normalize to [0, 1]
    surface -= surface.min()
    surface /= surface.max()

    return surface


# ============================================================
# Recurrence entropy
# ============================================================

def calculate_recurrence_entropy(image, threshold):
    """
    Calculate recurrence entropy and mean recurrence rate
    from a two-dimensional spatial pattern.

    Parameters
    ----------
    image : ndarray
        Two-dimensional image or spatial pattern.
    threshold : float
        Recurrence threshold.

    Returns
    -------
    entropy : float
        Normalized recurrence entropy.
    rho : float
        Mean recurrence rate.
    """

    # Moore neighborhood
    neighborhood = np.array([
        [-1, -1],
        [-1,  0],
        [-1,  1],
        [ 0,  1],
        [ 1,  1],
        [ 1,  0],
        [ 1, -1],
        [ 0, -1]
    ])

    n_neighbors = len(neighborhood)
    n_states = 2**n_neighbors

    # Powers used to encode the recurrence microstates
    powers = 2**np.arange(n_neighbors)

    # Statistics of recurrence microstates
    stats = np.zeros(n_states, dtype=np.int64)

    # Mean local recurrence rate
    rho_array = np.zeros_like(image)

    n_rows, n_cols = image.shape

    # Avoid boundary pixels
    for i in range(1, n_rows - 1):
        for j in range(1, n_cols - 1):

            center = image[i, j]
            state = 0
            recurrence_count = 0

            for k, (di, dj) in enumerate(neighborhood):

                neighbor = image[i + di, j + dj]

                recurrent = int(abs(center - neighbor) <= threshold)

                state += recurrent * powers[k]
                recurrence_count += recurrent

            stats[state] += 1

            rho_array[i, j] = recurrence_count / n_neighbors

    # Probability distribution of recurrence microstates
    probabilities = stats / stats.sum()

    # Shannon entropy
    entropy = -np.sum(
        probabilities[probabilities > 0]
        * np.log(probabilities[probabilities > 0])
    )

    # Normalize entropy
    entropy /= np.log(n_states)

    # Mean recurrence rate
    rho = np.mean(rho_array[1:-1, 1:-1])

    return entropy, rho


# ============================================================
# Save results
# ============================================================

def save_results(filename, results):
    """
    Save recurrence entropy results to a .dat file.
    """

    header = (
        "H\tseed\tN\tthreshold\tentropy\trecurrence_rate"
    )

    np.savetxt(
        filename,
        results,
        fmt="%.6f\t%d\t%d\t%.6f\t%.8f\t%.8f",
        header=header,
        comments=""
    )


# ============================================================
# Plot
# ============================================================

def plot_results(
    image,
    H,
    seed,
    thresholds,
    entropies,
    recurrence_rates,
    filename
):
    """
    Generate a figure showing the spatial pattern,
    recurrence entropy, and recurrence rate.
    """

    fig, axes = plt.subplots(
        1, 3,
        figsize=(15, 5)
    )

    # --------------------------------------------------------
    # Spatial pattern
    # --------------------------------------------------------

    axes[0].imshow(
        image,
        cmap="gray",
        vmin=0,
        vmax=1
    )

    axes[0].set_title(
        rf"Fractional Brownian surface ($H={H:.2f}$)"
    )

    axes[0].axis("off")

    # --------------------------------------------------------
    # Recurrence entropy
    # --------------------------------------------------------

    axes[1].plot(
        thresholds,
        entropies,
        linewidth=2
    )

    axes[1].set_xlabel("Threshold")
    axes[1].set_ylabel("Recurrence entropy $S$")
    axes[1].set_title("Recurrence entropy")

    axes[1].grid(alpha=0.3)

    # --------------------------------------------------------
    # Recurrence rate
    # --------------------------------------------------------

    axes[2].plot(
        thresholds,
        recurrence_rates,
        linewidth=2
    )

    axes[2].set_xlabel("Threshold")
    axes[2].set_ylabel(r"Mean recurrence rate $\rho$")
    axes[2].set_title("Recurrence rate")

    axes[2].grid(alpha=0.3)

    # --------------------------------------------------------
    # Figure layout
    # --------------------------------------------------------

    fig.suptitle(
        rf"$H={H:.2f}$, seed={seed}",
        fontsize=14
    )

    fig.tight_layout()

    fig.savefig(
        filename,
        dpi=200,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# Main analysis
# ============================================================

def main():

    # --------------------------------------------------------
    # Parameters
    # --------------------------------------------------------

    N = 256
    seed = 42

    H_values = np.linspace(0, 1, 11)
    threshold_values = np.linspace(0, 1, 51)

    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    data_file = os.path.join(
        output_dir,
        "recurrence_entropy_results.dat"
    )

    # Store all results
    all_results = []

    # --------------------------------------------------------
    # Loop over Hurst exponents
    # --------------------------------------------------------

    for H in H_values:

        print(f"Processing H = {H:.2f}")

        # Generate spatial pattern
        image = generate_fbm_surface(
            N=N,
            H=H,
            seed=seed
        )

        entropies = []
        recurrence_rates = []

        # ----------------------------------------------------
        # Loop over thresholds
        # ----------------------------------------------------

        for threshold in threshold_values:

            entropy, rho = calculate_recurrence_entropy(
                image,
                threshold
            )

            entropies.append(entropy)
            recurrence_rates.append(rho)

            all_results.append([
                H,
                seed,
                N,
                threshold,
                entropy,
                rho
            ])

        # ----------------------------------------------------
        # Generate figure
        # ----------------------------------------------------

        figure_file = os.path.join(
            output_dir,
            f"results_H_{H:.2f}_seed_{seed:04d}.png"
        )

        plot_results(
            image=image,
            H=H,
            seed=seed,
            thresholds=threshold_values,
            entropies=entropies,
            recurrence_rates=recurrence_rates,
            filename=figure_file
        )

    # --------------------------------------------------------
    # Save all numerical results
    # --------------------------------------------------------

    all_results = np.array(all_results)

    save_results(
        data_file,
        all_results
    )

    print("\nAnalysis completed.")
    print(f"Results saved to: {data_file}")
    print(f"Figures saved to: {output_dir}/")


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":
    main()
