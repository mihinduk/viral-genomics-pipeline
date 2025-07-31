# After args parsing, add:
    # Determine sample name
    if args.sample_name:
        sample_name = args.sample_name
    else:
        # Extract from output prefix (e.g., SMS_5_DENV1_results/SMS_5_DENV1 -> SMS_5)
        base_name = os.path.basename(args.output_prefix)
        # Try to extract sample name before underscore and virus name
        parts = base_name.split('_')
        if len(parts) >= 2:
            sample_name = '_'.join(parts[:2])  # e.g., SMS_5
        else:
            sample_name = parts[0]
    
    print(f"Using sample name: {sample_name}")

# Update the SeqRecord creations:
# For consensus:
            id=f"{sample_name}_filtered_consensus",

# For proteins:
            id=f"{sample_name}_{gene}",
