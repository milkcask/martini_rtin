import config from "eslint-config-mourner";

export default [
    ...config,
    {
        rules: {
            'no-use-before-define': 0
        }
    }
];
