eye_snps = {
    "rs12913832", "rs1800407", "rs1129038",
    "rs12203592", "rs16891982", "rs1426654",
    "rs11636232", "rs12896399", "rs1393350",
    "rs4778138"
}

def parse_dna_file(file_path):
    snp_data = {}

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            parts = line.strip().split()
            if len(parts) < 4:
                continue

            rsid = parts[0]
            genotype = parts[3]

            if rsid in eye_snps:
                snp_data[rsid] = genotype

    return snp_data


def predict_eye_color(snp_data):
    eye_color_probabilities = {
        "blue": 0.0,
        "brown": 0.0,
        "green": 0.0,
        "hazel": 0.0
    }

    if snp_data.get("rs12913832") == "AA":
        eye_color_probabilities["blue"] += 0.7
        eye_color_probabilities["brown"] -= 0.5
    elif snp_data.get("rs12913832") == "GG":
        eye_color_probabilities["brown"] += 0.7
        eye_color_probabilities["blue"] -= 0.5

    if snp_data.get("rs1800407") == "AA":
        eye_color_probabilities["blue"] += 0.6
        eye_color_probabilities["brown"] -= 0.4
    elif snp_data.get("rs1800407") == "GG":
        eye_color_probabilities["brown"] += 0.6
        eye_color_probabilities["blue"] -= 0.4

    total_prob = sum(eye_color_probabilities.values())
    for color in eye_color_probabilities:
        eye_color_probabilities[color] /= total_prob

    predicted_eye_color = max(eye_color_probabilities, key=eye_color_probabilities.get)
    return predicted_eye_color, eye_color_probabilities
